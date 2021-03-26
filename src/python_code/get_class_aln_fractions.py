from Bio import SeqIO
import csv
import os
# get the average class percents per column in every alignment. 

def get_frequencies(length, alignment):
  freq_dict = [dict(A = 0, R = 0, N = 0, D = 0, C = 0, E = 0, Q = 0, G = 0, H = 0, I = 0, L = 0, K = 0, M = 0, F = 0, P = 0, S = 0, T = 0, W = 0, Y = 0, V = 0) for x in range(length)]

  sum_list = []
  for i in range(length):
    sum = 0
    for seq in alignment:
      if seq[i] != "-" and seq[i] != "X" and seq[i] != "Z" and seq[i] != "B":
        freq_dict[i][seq[i]] += 1
        sum += 1
    sum_list.append(sum)

  for col, sum in zip(freq_dict, sum_list):
    for aa in col:
      if sum != 0:
        col[aa] = col[aa]/sum
      else:
        col[aa] = 0

  return freq_dict

def get_class_freq(freq, length):
  
  aliphatic = ["G", "A", "V", "L", "M", "I"]
  polar = ["S", "T", "C", "N", "Q"]
  positive = ["K", "R", "H"]
  negative = ["D", "E"]
  aromatic = ["F", "Y", "W"]
  proline = "P"

  class_freq = [dict(aliphatic=0, polar=0, positive=0, negative=0, aromatic=0, proline=0) for x in range(length)]
  for i, col in enumerate(freq): # there are "length" columns
    for aa in col:
      if aa in aliphatic:
        class_freq[i]["aliphatic"] += col[aa]
      if aa in polar:
        class_freq[i]["polar"] += col[aa] 
      if aa in positive:
        class_freq[i]["positive"] += col[aa]
      if aa in negative:
        class_freq[i]["negative"] += col[aa]
      if aa in aromatic:
        class_freq[i]["aromatic"] += col[aa]
      if aa == proline:
        class_freq[i]["proline"] += col[aa]

  return class_freq

#---------------main-------------------------
input_path = "../../data/PSICOV/aln_fasta/"
output_path = "../../output/output_PSICOV/" 

with open(output_path + "class_percents.csv", 'w', newline='\n', encoding='utf-8') as csv_file:
  writer = csv.writer(csv_file)
  writer.writerow(['gene', 'position', 'aa_class', 'frac_of_col'])

  classes = ["aliphatic", 
  protein_list = os.listdir(input_path)

  for protein in protein_list:
    records = list(SeqIO.parse(input_path + protein, "fasta"))
    ref_seq = records[0].seq

    
    for i in range(1, len(records)):
      aln_seq = records[i].seq
      percent_sum += get_percent(str(aln_seq), str(ref_seq))

    percent = percent_sum/len(records)
    writer.writerow([protein[0:4], percent])