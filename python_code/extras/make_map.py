#!/usr/bin/python

'''
This script parses an input PDB file and FASTA sequence and maps the sequence
to the PDB file using mafft. The script produces a CSV with columns for with
PDB residue numbering, corresponding FASTA sequence numbering, and the amino
acid. The map produced will have gaps because not all amino acids may be
respresented in teh PDB structure.

Author: Benjamin R. Jack
'''

import os
import csv
import argparse
import tempfile
import textwrap
import subprocess

from Bio import SeqIO, AlignIO
from Bio.Data import SCOPData
from Bio.PDB import PDBParser, is_aa
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_aa_seq(chain):
    '''
    Extract amino acid sequence from a PDB chain object and return sequence as
    Bio.SeqRecord object.
    '''
    aa_list = []
    residue_numbers = []
    for residue in chain:
        if is_aa(residue):
            aa_list.append(SCOPData.protein_letters_3to1[residue.resname])
            residue_numbers.append(str(residue.get_id()[1]) +
                                   residue.get_id()[2].strip())
    aa_seq = SeqRecord(Seq(''.join(aa_list)), id='pdb_seq', description='')
    return aa_seq, residue_numbers


def run_mafft(fasta_aln, pdb_seq):
    '''
    Align two Bio.SeqRecord sequences with mafft and return an a biopython 
    alignment object.
    '''
    # Open temporary files for mafft
    temp_fasta_pdb = tempfile.NamedTemporaryFile()
    temp_aln_clipped = tempfile.NamedTemporaryFile()
    temp_aln_out = tempfile.NamedTemporaryFile()
    # Temporary file for mafft output
    SeqIO.write(pdb_seq, temp_fasta_pdb.name, "fasta")
    try:
        print("Running mafft...\n")
        # Align sequences while keeping length
        subprocess.call(['mafft-linsi', '--maxiterate', '0',
                         '--keeplength', '--add',
                         temp_fasta_pdb.name, fasta_aln],
                        stdout=temp_aln_clipped)
        # Realign to see if pdb seq has been clipped on the ends
        subprocess.call(['mafft-linsi', '--add',
                         temp_fasta_pdb.name, temp_aln_clipped.name],
                        stdout=temp_aln_out)
    except:
        raise RuntimeError('Call to mafft failed. Check that mafft is '
                           'in your PATH.')
    alignment = AlignIO.read(temp_aln_out.name, 'fasta')
    temp_fasta_pdb.close()
    temp_aln_clipped.close()
    temp_aln_out.close()
    return alignment


def load_pdb_chain(name, pdb_file, model_name, chain_name):
    '''
    Load a specified chain from a PDB, with error checking.
    '''
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(name, pdb_file)
    try:
        model = structure[model_name]
    except KeyError:
        raise RuntimeError('PDB model could not be found. Please inspect PDB'
                           'file and specify model for map.')
    try:
        chain = model[chain_name]
    except KeyError:
        raise RuntimeError('PDB chain could not be found. Please inspect PDB'
                           'file and specify chain for map.')
    return chain


def make_map(alignment, residue_numbers, chain_name):
    '''
    Make a amino acid to PDB residue number map using an alignment and a list of
    PDB residue numbers. Return a list of dictionaries designed to be converted
    to a CSV.
    '''
    # Split aligned sequences into two lists
    pdb_clipped = list(alignment[-2])
    pdb_full = list(alignment[-1])
    # Track *index* of where we are in the PDB amino acid sequence
    pdb_index = 0
    # Track alignment position (starts at 1, not an index!)
    aln_position = 1
    out_list = []
    for pdb_clipped_aa, pdb_full_aa in zip(pdb_clipped, pdb_full):
        out_dict = {}
        if pdb_clipped_aa == '-' and pdb_full_aa == '-':
            # This is a gap in the original alignment
            out_dict['pdb_position'] = 'NA'
            out_dict['pdb_aa'] = 'NA'
            out_dict['chain'] = 'NA'
            out_dict['aln_position'] = aln_position
            aln_position += 1
        elif pdb_clipped_aa == '-' and pdb_full_aa != '-':
            # Part of pdb sequence must have been clipped off
            out_dict['pdb_position'] = residue_numbers[pdb_index]
            out_dict['pdb_aa'] = pdb_full_aa
            out_dict['chain'] = chain_name
            pdb_index += 1
            out_dict['aln_position'] = 'NA'
        elif pdb_clipped_aa != '-' and pdb_full_aa != '-':
            out_dict['pdb_position'] = residue_numbers[pdb_index]
            out_dict['pdb_aa'] = pdb_full_aa
            out_dict['chain'] = chain_name
            pdb_index += 1
            out_dict['aln_position'] = aln_position
            aln_position += 1
        else:
            raise RuntimeError(
                "The full PDB sequence contains a gap where the clipped PDB sequence does not. PDB sequence has mis-aligned to itself.")
        out_list.append(out_dict)

    return out_list


def main():
    '''
    Make a PDB amino acid to FASTA amino acid sequence map. This map is required
    to align evolutionary rates to positions in the structure.
    '''
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Generate sequence-to-structure map for aligning '
                    'evolutionary rates to a PDB structure.',
        epilog=textwrap.dedent('''\
            This script produces a CSV with the following columns, where each 
            row is a position in the sequence map:

            Column name     Description
            ===================================================================
            aln_position    Numeric position (column) in the multiple sequence
                            alignment extracted from FASTA file.
            
            pdb_position    Numeric position in sequence extracted from PDB 
                            file. This value corresponds to the numbered 
                            positions generated by `calc_wcn.py` and 
                            `calc_rsa.py`.

            pdb_aa          Single letter amino acid for this position extracted
                            from the PDB file. Missing values are positions not 
                            present in the original PDB file.

            chain           PDB chain corresponding to this position.

            ''')
    )
    parser.add_argument('fasta', metavar='<FASTA alignment>', type=str,
                        help='input FASTA multiple sequence alignment')
    parser.add_argument('pdb', metavar='<PDB path>', type=str,
                        help='input PDB file')
    parser.add_argument('-c', metavar='<PDB chain>', type=str,
                        default='A',
                        help='if there are multiple chains in PDB, map FASTA '
                             'sequence to this chain')
    parser.add_argument('-m', metavar='<PDB model>', type=int,
                        default=0,
                        help='if there are multiple models in PDB, map FASTA '
                             'sequence to this model')
    parser.add_argument('-o', metavar='<output prefix>', type=str,
                        help='prefix for output files')
    args = parser.parse_args()
    # Grab PDB name from filename
    pdb_name = os.path.splitext(os.path.basename(args.pdb))[0]
    # Define output file names
    if args.o is None:
        # If no output prefix given, assign prefix using input filename
        args.o = pdb_name
    output_map = args.o + '.map.csv'
    # Load chain
    print("Extracted amino acid sequence from chain '{}' in model '{}' of PDB "
          "file '{}'.".format(args.c, args.m, args.pdb))
    print("If you would like to specify a different chain or model, please "
          "run `make_map.py -h` for help.\n")
    chain = load_pdb_chain(pdb_name.upper(), args.pdb, args.m, args.c)
    # Extract PDB numbering and amino acid sequence
    pdb_record, residue_numbers = get_aa_seq(chain)
    # Align PDB to FASTA alignment
    alignment = run_mafft(args.fasta, pdb_record)
    print(alignment)
    # Generate map
    output_list = make_map(alignment, residue_numbers, args.c)
    # Write map to CSV
    with open(output_map, 'w') as csvfile:
        writer = csv.DictWriter(csvfile,
                                fieldnames=['aln_position',
                                            'pdb_position', 'pdb_aa', 'chain'],
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output_list)
    print("\nMap successfully generated.")


if __name__ == "__main__":
    main()
