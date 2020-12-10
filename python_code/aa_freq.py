import csv
import math
import os
import sys

#This program takes in multiple protein alignments and calculates the entropy for each aa position
'''
1) calc aa freq (20 values) for each col in alignment
2) cal entropy for each col (H). 
3) calc e^H. record in txt for each col/site
4) repeat for 38 proteins (aligned)
'''

#list of protein alignments
protein_list = os.listdir(path="./aligned_sequences/")

# output file with all n-eff values for each of 38 proteins. 
#result = open("variability.txt","w+")

with open("natural_variability_2.csv", "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file) 
  writer.writerow(['position', 'gene', 'q_H', 'q_E', 'q_D',  'q_R', 'q_K', 'q_S', 'q_T', 'q_N', 'q_Q', 'q_A', 'q_V', 'q_L', 'q_I', 'q_M', 'q_F', 'q_Y', 'q_W', 'q_P', 'q_G', 'q_C', 'entropy', 'n_eff'])

  for protein in protein_list:
    alignment = []
    with open("./aligned_sequences/" + protein, 'r') as file:
      lines = [line.rstrip('\n') for line in file]

    seq = ""
    for line in lines:
      if ">" in line and len(seq) != 0:
        alignment.append(seq)
        seq = ""
      elif ">" not in line:
        seq = seq + line
    if len(seq) != 0:
      alignment.append(seq)
    
    name = lines[0]
    #result.write(name + "\n")

    #length of the sequence (columns)
    length = len(alignment[0])
    sum_list = []

    # calulating freq of each aa by position in protein seq (2D array)
    freq = [dict(A = 0, R = 0, N = 0, D = 0, C = 0, E = 0, Q = 0, G = 0, H = 0, I = 0, L = 0, K = 0, M = 0, F = 0, P = 0, S = 0, T = 0, Y = 0, W = 0, V = 0) for x in range(length)]

    for i in range(length):
      s = 0
      for seq in alignment:
        if seq[i] != "-" and seq[i] != "X" and seq[i] != "Z" and seq[i] != "B":
          freq[i][seq[i]] += 1
          s += 1
      sum_list.append(s)

    for col, s in zip(freq, sum_list):
      for aa in col:
        col[aa] = col[aa]/s

    # calculating entropy and n-eff for each site
    n_eff = []
    entropy = [] 

    for col in freq:
      sum = 0
      for aa in col:
        sum += col[aa] * (1 if col[aa] == 0 else math.log(col[aa]))
      entropy.append(-sum)

    for value in entropy:
      n_eff.append(math.exp(value))

    #result.write(str(n_eff) + "\n")

    aa = ['H', 'E', 'D',  'R', 'K', 'S', 'T', 'N', 'Q', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'P', 'G', 'C']
    #position = 0
    for position in range(length):
      freq_list = []
      for key in aa:
        freq_list.append(freq[position][key])
      writer.writerow([str(position + 1), str(protein[0:4]).lower()] + freq_list + [str(entropy[position]), str(n_eff[position])]) 

#file.close()


