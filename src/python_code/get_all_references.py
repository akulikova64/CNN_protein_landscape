from Bio import SeqIO
import os
import sys


# get the reference into fasta

input_path = "../../data/PSICOV/aln_fasta/"
output_path = "../../data/PSICOV/PSICOV_seqs.fasta"

protein_list = os.listdir(input_path)

with open(output_path, 'w') as file:
  for protein in protein_list:
    records = list(SeqIO.parse(input_path + protein, "fasta"))
    ref_seq = records[0].seq

    file.write(">" + str(protein[0:4]) + "\n")
    file.write(str(ref_seq) + "\n")