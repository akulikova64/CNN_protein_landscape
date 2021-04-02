from Bio import SeqIO
import csv
import sys
import os

# get cleam fasta format

input = "../../data/PSICOV/new_training_seqs.fasta" # csv with all stats
output = "../../data/PSICOV/new_training_seqs_clean.fasta" # csv with all stats

with open(output, 'w', newline='\n', encoding='utf-8') as file:
 
  records = list(SeqIO.parse(input, "fasta")) #PSICOV record object
  for record in records:
    file.write(">" + str(record.name) + "\n")
    file.write(str(record.seq) + "\n")