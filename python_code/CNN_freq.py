import csv
import math
import sys
import os

# This program takes in aa acid frequencies at each position from raw CNN output. 

'''
1) Calculate the entropy and n-eff for each position
2) Save to csv file
'''

def get_entropy(freq_list):
  sum = 0
  for freq in freq_list:
    freq = float(freq)
    sum += freq * (1 if freq == 0 else math.log(freq))
  return -sum

def get_n_eff(entropy):

  return math.exp(entropy)


# --------------- main ------------------------------------
output_path = "../data/output/stats_cnn.csv"
input_path = "../data/CNN_output/"


gene_files = os.listdir(path="./CNN_output_new")

# output CSV file with all the CNN entropy, n-eff and freq values for each of 38 proteins. 
with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file) 
  writer.writerow(['position', 'gene', 'wt_aa', 'q_H', 'q_E', 'q_D',  'q_R', 'q_K', 'q_S', 'q_T', 'q_N', 'q_Q', 'q_A', 'q_V', 'q_L', 'q_I', 'q_M', 'q_F', 'q_Y', 'q_W', 'q_P', 'q_G', 'q_C', 'entropy', 'n_eff'])

  # freq of aa acids (changes for each position)
  for gene in gene_files:
    with open(input_path + gene, 'r') as file:
      gene = gene[:-4]
      file.readline()
      position = 0
      for line in file:
        line = line.rstrip('\n')
        line = line.rsplit()
        position += 1
        wt_aa = line[1]
        
        aa = ['H', 'E', 'D',  'R', 'K', 'S', 'T', 'N', 'Q', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'P', 'G', 'C']
        freq = {'H':0, 'E':0, 'D':0, 'R':0, 'K':0, 'S':0, 'T':0, 'N':0, 'Q':0, 'A':0, 'V':0, 'L':0, 'I':0, 'M':0, 'F':0, 'Y':0, 'W':0, 'P':0, 'G':0, 'C':0}
       
        for i in range(2, 22):
          freq[aa[i-2]] = line[i]

        # get 20 frequencies from "line"
        freq_list = []
        for key in aa:
          freq_list.append(freq[key])

        entropy = get_entropy(freq_list)
        n_eff = get_n_eff(entropy)

        writer.writerow([str(position), str(gene[0:4]).lower(), str(wt_aa)] + freq_list + [str(entropy), str(n_eff)]) 

