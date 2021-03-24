import csv
import math
import os
import sys

#This program takes in multiple sequence alignments and calculates the entropy for each aa position
'''
1) Calculates aa freq (20 values) for each col in alignment
2) Calculates entropy for each col (H)
3) Calculates e^H
4) Records into CSV for each col/site
4) Repeats for all proteins (aligned)
'''

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

def get_frequencies(length, alignment):
  freq = [dict(A = 0, R = 0, N = 0, D = 0, C = 0, E = 0, Q = 0, G = 0, H = 0, I = 0, L = 0, K = 0, M = 0, F = 0, P = 0, S = 0, T = 0, W = 0, Y = 0, V = 0) for x in range(length)]

  sum_list = []
  for i in range(length):
    sum = 0
    for seq in alignment:
      if seq[i] != "-" and seq[i] != "X" and seq[i] != "Z" and seq[i] != "B":
        freq[i][seq[i]] += 1
        sum += 1
    sum_list.append(sum)

  for col, sum in zip(freq, sum_list):
    for aa in col:
      if sum != 0:
        col[aa] = col[aa]/sum
      else:
        col[aa] = 0

  return freq

def get_entropy(freq):
  entropy = []
  for col in freq:
      sum = 0
      for item in col: # item can be aa or aa_class
        sum += col[item] * (1 if col[item] == 0 else math.log(col[item]))
      entropy.append(-sum)

  return entropy

def get_n_eff(entropy):

  n_eff = []
  for value in entropy:
      n_eff.append(math.exp(value))
  
  return n_eff

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

# --------- main -------------------------
group_name = "100"
min_seqs = "10" #minimum sequences per alignment (choose out of the ones made)
output_path = "../../output/output_PSICOV/stats_align_files/stats_align_" + group_name + ".csv"

# for group names: 40, 60 and 80, alignments must be longer than "min_seqs". 
if group_name != "20" or group_name != "100":
  input_path = "../../data/PSICOV/aln_filtered_" + str(min_seqs) + "/aln_" + group_name + "/" #input multiple sequence aligments path
else:
  input_path = "../../data/PSICOV/aln_" + group_name + "/" #input multiple sequence aligments path

#list of protein sequence alignments:
protein_list = os.listdir(input_path)

with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file) 
  # old: writer.writerow(['position', 'gene', 'q_H', 'q_E', 'q_D',  'q_R', 'q_K', 'q_S', 'q_T', 'q_N', 'q_Q', 'q_A', 'q_V', 'q_L', 'q_I', 'q_M', 'q_F', 'q_Y', 'q_W', 'q_P', 'q_G', 'q_C', 'entropy', 'n_eff', 'q_aliphatic', 'q_polar', 'q_positive', 'q_negative', 'q_aromatic', 'q_proline', 'entropy_class', 'n_eff_class'])
  writer.writerow(['position', 'gene', 'q_A', 'q_R', 'q_N',  'q_D', 'q_C', 'q_Q', 'q_E', 'q_G', 'q_H', 'q_I', 'q_L', 'q_K', 'q_M', 'q_F', 'q_P', 'q_S', 'q_T', 'q_W', 'q_Y', 'q_V', 'entropy', 'n_eff', 'q_aliphatic', 'q_polar', 'q_positive', 'q_negative', 'q_aromatic', 'q_proline', 'entropy_class', 'n_eff_class'])

  for protein in protein_list:
    with open(input_path + protein, 'r') as file:
      lines = [line.rstrip('\n') for line in file]

    if len(lines) == 0:
      continue

    # get a list of aligned sequences
    alignment = get_alignment(lines)
    name = lines[0]

    # get length of the sequence (# of columns):
    length = len(alignment[0])

    # calulate freq of each aa by position in protein seq (2D array):
    freq = get_frequencies(length, alignment)

    # calculate entropy and n-eff for each site:
    entropy = get_entropy(freq)
    n_eff = get_n_eff(entropy)

    # get 6 class frequencies from class_freq (for each site)
    class_freq = get_class_freq(freq, length)

    entropy_class = get_entropy(class_freq)
    n_eff_class = get_n_eff(entropy_class)

    # saving freq, entropy and n_eff values to CSV:
    #old: aa = ['H', 'E', 'D', 'R', 'K', 'S', 'T', 'N', 'Q', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'P', 'G', 'C']
    aa = ['A', 'R', 'N',  'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    classes = ['aliphatic', 'polar', 'positive', 'negative', 'aromatic', 'proline']
    
    for position in range(length):
      freq_list = []
      class_freq_list = []
      for key in aa:
        freq_list.append(freq[position][key])
      for key in classes:
        class_freq_list.append(class_freq[position][key])

      writer.writerow([str(position + 1), str(protein[0:4]).lower()] + freq_list + [str(entropy[position]), str(n_eff[position])] + class_freq_list + [str(entropy_class[position]), str(n_eff_class[position])]) 

print("Saved CSV to " + output_path)



