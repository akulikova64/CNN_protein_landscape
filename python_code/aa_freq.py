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
4) Repeats for 38 proteins (aligned)
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
  freq = [dict(A = 0, R = 0, N = 0, D = 0, C = 0, E = 0, Q = 0, G = 0, H = 0, I = 0, L = 0, K = 0, M = 0, F = 0, P = 0, S = 0, T = 0, Y = 0, W = 0, V = 0) for x in range(length)]

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
      col[aa] = col[aa]/sum

  return freq

def get_entropy(freq):
  entropy = []
  for col in freq:
      sum = 0
      for aa in col:
        sum += col[aa] * (1 if col[aa] == 0 else math.log(col[aa]))
      entropy.append(-sum)

  return entropy

def get_n_eff(entropy):

  n_eff = []
  for value in entropy:
      n_eff.append(math.exp(value))
  
  return n_eff


# --------- main -------------------------
output_path = "../data/output/natural_var_align_all.csv"
input_path = "../data/data_all/"


#list of protein sequence alignments:
protein_list = os.listdir(input_path)

with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file) 
  writer.writerow(['position', 'gene', 'q_H', 'q_E', 'q_D',  'q_R', 'q_K', 'q_S', 'q_T', 'q_N', 'q_Q', 'q_A', 'q_V', 'q_L', 'q_I', 'q_M', 'q_F', 'q_Y', 'q_W', 'q_P', 'q_G', 'q_C', 'entropy', 'n_eff'])

  for protein in protein_list:
    with open(input_path + protein, 'r') as file:
      lines = [line.rstrip('\n') for line in file]

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

    # saving freq, entropy and n_eff values to CSV:
    aa = ['H', 'E', 'D', 'R', 'K', 'S', 'T', 'N', 'Q', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'P', 'G', 'C']
    #position = 0
    for position in range(length):
      freq_list = []
      for key in aa:
        freq_list.append(freq[position][key])
      writer.writerow([str(position + 1), str(protein[0:4]).lower()] + freq_list + [str(entropy[position]), str(n_eff[position])]) 




