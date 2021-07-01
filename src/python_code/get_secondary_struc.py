import csv
import math
import os
import sys

# extracting the secondary structure from the STRIDE output files. 

input_path = "../../data/PSICOV/secondary_structure/" #input multiple sequence aligments path
output_path = "../../output/output_PSICOV/second_struc.csv"

#list of STRIDE files:
protein_list = os.listdir(input_path)


with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file) 
  # old: writer.writerow(['position', 'gene', 'q_H', 'q_E', 'q_D',  'q_R', 'q_K', 'q_S', 'q_T', 'q_N', 'q_Q', 'q_A', 'q_V', 'q_L', 'q_I', 'q_M', 'q_F', 'q_Y', 'q_W', 'q_P', 'q_G', 'q_C', 'entropy', 'n_eff', 'q_aliphatic', 'q_polar', 'q_positive', 'q_negative', 'q_aromatic', 'q_proline', 'entropy_class', 'n_eff_class'])
  writer.writerow(['gene', 'position', 'second_struc'])

  for protein in protein_list:
    with open(input_path + protein, 'r') as file:
      gene = protein[0:4]
      for line in file:
        line = line.strip('\n')
        line = line.split()

        if line[0] == "ASG":
          position = line[3]
          sec_struc = line[6]
          writer.writerow([str(gene), int(position), str(sec_struc)]) 

print("Saved CSV to " + output_path)