from Bio import SeqIO
import os
import sys

# This code samples the sequences PSICOV alignments by % similarity to the sequences in the PSICOV pdb structures.

#Plan:
#load the seqs from the CNN output (wt aa) into a dictionary with the protein name. 
#open and parse each of  the fasta files 
# make sure that name of fasta file matches with name the pdb sequences. 

#loading a fasta aln file:

MAX_LIM = 0.4
MIN_LIM = MAX_LIM - 0.2

def compare(aln_seq, ref_seq):
  count_similar = 0
  for i, j in zip(aln_seq, ref_seq):
    assert j != "-"
    if i == j:
      count_similar += 1

  similarity = float(count_similar/len(ref_seq))

  if similarity > MIN_LIM and similarity <= MAX_LIM:
    return True
  else:
    return False

# do not use this function (it is just in case)
def get_reference_from_CNN_data():
  
  aaCodes = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T','TRP':'W', 'TYR':'Y', 'VAL':'V'}
  with open("../data/PSICOV/PSICOV_CNN_output/1a6m_final_tot.csv", 'r') as file:
    ref_seq = ""
    file.readline()
    for line in file:
      line = line.rstrip('\n')
      line = line.rsplit(",")
      wt_aa = line[5]
      ref_seq += aaCodes[wt_aa]

  return ref_seq

#CNN_reference = get_reference_from_CNN_data() #should be the same as records[0] in alignment

input_path = "../../data/PSICOV/aln_fasta/"
output_path = "../../data/PSICOV/" + "aln_" + str(int(MAX_LIM*100)) + "/"

protein_list = os.listdir(input_path)

for protein in protein_list:
  with open(output_path + protein, 'w') as file:
    records = list(SeqIO.parse(input_path + protein, "fasta"))
    ref_seq = records[0].seq
    #file.write(">reference\n") # we do not want to include the reference in the calculations
    #file.write(str(ref_seq) + "\n")

    for i in range(1, len(records)):
      aln_seq = records[i].seq
      keep = compare(str(aln_seq), str(ref_seq))

      if keep:
        file.write(">" + str(i) + "\n")
        file.write(str(aln_seq) + "\n")

      

