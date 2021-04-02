from Levenshtein import distance as levenshtein_distance
from Bio import SeqIO
import csv
import sys
import os

def calc_sim_corrected(seq1, seq2):

    min_len = min(len(seq1), len(seq2))
    abs_diff = abs(len(seq1) - len(seq2))
    lev_dist  = levenshtein_distance(seq1, seq2) - abs_diff
    return (1 - lev_dist / min_len)*100, lev_dist

# calculating the levenshtein distance between sequence pairs:

# fasta files:
'''
fasta_1 = "../../data/PSICOV/PSICOV_seqs.fasta" # analysis sequences
fasta_2 = "../../data/PSICOV/new_training_seqs.fasta" # CNN training sequeces'''

# testing sequences:

fasta_1 = "./test_seq1.fasta"
fasta_2 = "./test_seq2.fasta"

output_path_1 = "../../data/PSICOV/above_50_similarity_test.csv" # csv with all stats
output_path_2 = "../../data/PSICOV/above_50_similarity_test.txt" # only the pdb ids and chains

# creating a csv with all similarity and sequence stats:
with open(output_path_1, 'w', newline='\n', encoding='utf-8') as csv_file:
  writer = csv.writer(csv_file) 
  writer.writerow(['CNN_seq_id', "CNN_seq_length", "absolute_difference", "lev_dist_original", "perc_sim_original", "lev_dist_corrected", "perc_sim_corrected", "PSICOV_seq_id", "PSICOV_seq_length"])

  records_1 = list(SeqIO.parse(fasta_1, "fasta")) #PSICOV record object
  records_2 = list(SeqIO.parse(fasta_2, "fasta")) # CNN training data object
  seq_list = []

  for rec_1 in records_1:
    for rec_2 in records_2:
    
      seq_1 = str(rec_1.seq) # PSICOV sequence
      seq_2 = str(rec_2.seq) # CNN training data sequence

      # calculating the levenshtein distance:
      lev_dist_original = levenshtein_distance(seq_1, seq_2) 

      name_1 = rec_1.name #PSICOV protien name
      name_2 = rec_2.name #CNN protein name

      # converting levenshtein distance to similarity
      similarity_corrected, lev_dist_corrected = calc_sim_corrected(seq_1, seq_2)
      
      similarity_original = (1 - (lev_dist_original/max(len(seq_1), len(seq_2))))*100
  
      # recording only the cases where the % similarity is greater than 50%:
      if similarity_original > 50 or similarity_corrected > 50:
        writer.writerow([str(name_2), str(len(seq_2)), str(abs(len(seq_1) - len(seq_2))), str(lev_dist_original), str(similarity_original), str(lev_dist_corrected), str(similarity_corrected), str(name_1), str(len(seq_1))])
        if name_2 not in seq_list:
          seq_list.append(name_2)

  # creating a txt file with only protein ID's (of protein that are greater than 50% similar)
  with open(output_path_2, 'w', newline='\n', encoding='utf-8') as txt_file:       
    for seq in seq_list:
      txt_file.write(str(seq) + "\n")  
      


