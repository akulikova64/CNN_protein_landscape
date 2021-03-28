from Levenshtein import distance as levenshtein_distance
from Bio import SeqIO
import csv
import sys
import os

# calculating the levenshtein distance between sequence pairs:

# fasta files:
fasta_1 = "../../data/PSICOV/PSICOV_seqs.fasta" # analysis sequences
fasta_2 = "../../data/PSICOV/training_seqs.fasta" # CNN training sequeces

'''
fasta_1 = "./test_seq1.fasta"
fasta_2 = "./test_seq2.fasta"'''

output_path = "../../data/PSICOV/above_50_similarity.csv"

with open(output_path, 'w', newline='\n', encoding='utf-8') as csv_file:
  writer = csv.writer(csv_file) 
  writer.writerow(['CNN_seq_id', "CNN_seq_length", "lev_dist", "perc_sim", "PSICOV_seq_id", "PSICOV_seq_length"])

  records_1 = list(SeqIO.parse(fasta_1, "fasta"))
  records_2 = list(SeqIO.parse(fasta_2, "fasta"))
  seq_list = []

  for rec_1 in records_1:
    for rec_2 in records_2:
    
      seq_1 = str(rec_1.seq)
      seq_2 = str(rec_2.seq)

      lev_dist = levenshtein_distance(seq_1, seq_2)

      name_1 = rec_1.name #PSICOV protien name
      name_2 = rec_2.name #CNN protein name

      similarity = (1 - lev_dist / max(len(seq_1), len(seq_2))) * 100

      if similarity > 50:
        writer.writerow([str(name_2), str(len(seq_2)), str(lev_dist), str(similarity), str(name_1), str(len(seq_1))])
        if name_2 not in seq_list:
          seq_list.append(name_2)
  '''    
  for seq in seq_list:
    file.write(str(seq) + "\n")  '''
      


