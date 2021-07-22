# import Bio import 
import freesasa
import get_freeSASA_scores
import os
import sys
import csv

# solvent_accessibility using freeSASA
aaCodes = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T','TRP':'W', 'TYR':'Y', 'VAL':'V'}

input_path = "../../../data/PSICOV/pdb/" #input multiple sequence aligments path
output_path = "../../../output/output_PSICOV/SASA_scores.csv"

# list of PDB files:
protein_list = os.listdir(input_path)

with open(output_path, "w+", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file) 
  writer.writerow(['gene', 'position', 'aa', 'SASA_abs_total', 'SASA_rel_total'])

  # parsing through all proteins
  for protein in protein_list:
    gene = protein[0:4]
    structure = freesasa.Structure(input_path + protein)
    result = freesasa.calc(structure)
    residues = result.residueAreas()

    # parsing though all residues
    for chain in residues:
      for key in residues[chain]:
        residue = residues[chain][key]
        position = residue.residueNumber
        aa = aaCodes[residue.residueType]
        sasa_abs = residue.total
        sasa_rel = residue.relativeTotal

        writer.writerow([str(gene), int(position), str(aa), str(sasa_abs), str(sasa_rel)]) 

print("Saved CSV to " + output_path)



