import csv
import math
import os
import sys

def get_alignment(lines):
  alignment = []
  seq = ""
  for line in lines:
    if ">" in line and len(seq) != 0:
      alignment.append(seq)
      seq = ""
    elif ">" not in line:
      seq = seq + line
  if len(seq) != 0:
    alignment.append(seq)

  return alignment


input_path = "../../data/PSICOV/aln_fasta/" #input multiple sequence aligments path
output_path = "../../output/output_PSICOV/natural_abundance.csv"

#list of protein sequence alignments:
protein_list = os.listdir(input_path)

with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file) 
  # old: writer.writerow(['position', 'gene', 'q_H', 'q_E', 'q_D',  'q_R', 'q_K', 'q_S', 'q_T', 'q_N', 'q_Q', 'q_A', 'q_V', 'q_L', 'q_I', 'q_M', 'q_F', 'q_Y', 'q_W', 'q_P', 'q_G', 'q_C', 'entropy', 'n_eff', 'q_aliphatic', 'q_polar', 'q_positive', 'q_negative', 'q_aromatic', 'q_proline', 'entropy_class', 'n_eff_class'])
  writer.writerow(['aa', 'count'])

  counts = {'A':0, 'R':0, 'N':0, 'D':0, 'C':0, 'E':0, 'Q':0, 'G':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0}

  for protein in protein_list:
    with open(input_path + protein, 'r') as file:
      lines = [line.rstrip('\n') for line in file]

    if len(lines) == 0:
      continue

    # get a list of aligned sequences
    alignment = get_alignment(lines)

    # get length of the sequence (# of columns):
    length = len(alignment[0])

    # calulate counts:
    for i in range(length):
      for seq in alignment:
        if seq[i] != "-" and seq[i] != "X" and seq[i] != "Z" and seq[i] != "B":
          counts[seq[i]] += 1

  # writing counts to CSV:
  for key, value in counts.items():
    aa = key
    count = value

    writer.writerow([str(aa), int(value)]) 

print("Saved CSV to " + output_path)