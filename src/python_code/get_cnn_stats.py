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

def get_class_freq(freq):
  
  unique = ["P", "G"]
  aliphatic = ["M", "L", "I", "V", "A"]
  small_polar = ["C", "S", "T", "N", "Q"]
  negative = ["D", "E"]
  positive = ["R", "K"]
  aromatic = ["H", "Y", "F", "W"]

  class_dict = {"unique":0, "aliphatic":0, "small polar":0, "negative":0, "positive":0, "aromatic":0}
  for aa in freq:
    if aa in unique:
      class_dict["unique"] += freq[aa]
    if aa in aliphatic:
      class_dict["aliphatic"] += freq[aa] 
    if aa in small_polar:
      class_dict["small polar"] += freq[aa]
    if aa in negative:
      class_dict["negative"] += freq[aa]
    if aa in positive:
      class_dict["positive"] += freq[aa]
    if aa == aromatic:
      class_dict["aromatic"] += freq[aa]

    class_freq_list = []
    for key in class_dict:
      class_freq_list.append(class_dict[key])

  return class_dict, class_freq_list

# --------------- main ------------------------------------
box_size_list = ["12", "20", "30", "40"]

for box_size in box_size_list:

  output_path = "../../data/PSICOV_box_" + box_size + "/output/stats_cnn.csv"
  input_path = "../../data/PSICOV_box_" + box_size + "/PSICOV_CNN_output/"


  gene_files = os.listdir(input_path)

  # output CSV file with all the CNN entropy, n-eff and freq values for each of 38 proteins. 
  with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:
    writer = csv.writer(CSV_file) 
    writer.writerow(['position', 'gene', 'wt_aa', 'q_A', 'q_R', 'q_N',  'q_D', 'q_C', 'q_Q', 'q_E', 'q_G', 'q_H', 'q_I', 'q_L', 'q_K', 'q_M', 'q_F', 'q_P', 'q_S', 'q_T', 'q_W', 'q_Y', 'q_V', 'entropy', 'n_eff', 'q_aliphatic', 'q_polar', 'q_positive', 'q_negative', 'q_aromatic', 'q_proline', 'entropy_class', 'n_eff_class'])

    # freq of aa acids (changes for each position)
    for gene in gene_files:
      with open(input_path + gene, 'r') as file:
        gene = gene[:-4]
        file.readline()
        for line in file:
          line = line.rstrip('\n')
          line = line.rsplit(",")
          position = str(int(line[0]) + 1)
          wt_aa = line[5]
          
          # prALA,prARG,prASN,prASP,prCYS,prGLN,prGLU,prGLY,prHIS,prILE,prLEU,prLYS,prMET,prPHE,prPRO,prSER,prTHR,prTRP,prTYR,prVAL
          aa = ['A', 'R', 'N',  'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
          freq = {'A':0, 'R':0, 'N':0, 'D':0, 'C':0, 'Q':0, 'E':0, 'G':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0}
        
          for i in range(10, 30):
            freq[aa[i-10]] = float(line[i])

          # get 20 frequencies from "line"
          freq_list = []
          for key in aa:
            freq_list.append(freq[key])
          
          entropy = get_entropy(freq_list)
          n_eff = get_n_eff(entropy)

          # get 6 class frequencies from class_freq
          class_freq, class_freq_list = get_class_freq(freq)

          entropy_class = get_entropy(class_freq_list)
          n_eff_class = get_n_eff(entropy_class)

          
          writer.writerow([str(position), str(gene[0:4]).lower(), str(wt_aa)] + freq_list + [str(entropy), str(n_eff)] + class_freq_list + [str(entropy_class), str(n_eff_class)]) 

  print("Saved CSV to " + output_path)

