import csv
import os

# count the number of seqs per gene alignment for each sequence similarity group:
groups = ["20","40","60","80","100"]
min_seqs = "10"

with open("../../output/output_PSICOV/seq_counts.csv", "w", newline='\n', encoding='utf-8') as counts_table:
  writer = csv.writer(counts_table) 
  writer.writerow(['group', 'gene', 'seq_count'])

  for group_name in groups:
    input_path = "../../data/PSICOV/aln_filtered_" + str(min_seqs) + "/aln_" + group_name + "/" #multiple sequence aligments path
    
    #list of protein sequence alignments:
    protein_list = os.listdir(input_path)

    for protein in protein_list:
      with open(input_path + protein, 'r') as file:
        lines = [line.rstrip('\n') for line in file]
        count = len(lines)/2

        writer.writerow([group_name, protein[0:4], count])

