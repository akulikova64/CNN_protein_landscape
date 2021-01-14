import os
import shutil
import csv
import sys

def getMax(list):
  aaList = ['H', 'E', 'D',  'R', 'K', 'S', 'T', 'N', 'Q', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'P', 'G', 'C']
  ind = 0
  max = 0
  for i in range(0,len(list)):
    if(float(list[i]) > max):
      ind = i
      max = float(list[i])

  return [aaList[ind], max]

def getMax_class(list):
  class_list = ["aliphatic", "polar", "positive", "negative", "aromatic", "proline"]
  ind = 0
  max = 0
  for i in range(0, len(list)):
    if(float(list[i]) > max):
      ind = i
      max = float(list[i])

  return [class_list[ind], max]

def get_aa_class(wt_aa, class_freqs):
  class_list = ["aliphatic", "polar", "positive", "negative", "aromatic", "proline"]

  aliphatic = ["G", "A", "V", "L", "M", "I"]
  polar = ["S", "T", "C", "N", "Q"]
  positive = ["K", "R", "H"]
  negative = ["D", "E"]
  aromatic = ["F", "Y", "W"]
  proline = "P"

  # get aa_class of the wt (from "wt_aa" parameter)  
  if wt_aa in aliphatic:
    aa_class = "aliphatic"
  if wt_aa in polar:
    aa_class = "polar" 
  if wt_aa in positive:
    aa_class = "positive"
  if wt_aa in negative:
    aa_class = "negative"
  if wt_aa in aromatic:
    aa_class = "aromatic"
  if wt_aa == proline:
    aa_class = "proline"

  ind = class_list.index(aa_class)
  class_freq = class_freqs[ind]

  return aa_class, class_freq

#--------- main ------------------------------------------
input_path_1 = "../../data/PSICOV/output/cnn_wt_max_freq.csv"
input_path_2 = "../../data/PSICOV/output/stats_align_all.csv"
output_path = "../../data/PSICOV/output/natural_max_freq.csv" 

wt = {}
# reading in CNN max frequencies file
with open(input_path_1, "r", newline='\n', encoding='utf-8') as CSV_file:
  csv_reader = csv.reader(CSV_file, delimiter=',')
  # gene, group, position, aa, freq
  #  0     1        2      3    4
  gene = ""
  for row in csv_reader:
    if row[1] == "wt": 
      if row[0] != gene: # making sure we do not create multiple dicts for the same gene.  
        wt[row[0]] = {} # starting a dict for the gene
        gene = row[0] # update current gene
      wt[row[0]][row[2]] = row[3] # wt[gene][position] = aa

with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file)
  writer.writerow(['gene', 'group', 'position', 'aa', 'freq', 'aa_class', 'class_freq'])

  # reading the alignment freq file
  with open(input_path_2, "r", newline='\n', encoding='utf-8') as CSV_file:
    csv_reader = csv.reader(CSV_file, delimiter=',')
    for row in csv_reader:
      if row[0] == 'position': # skipping header
        continue

      #*************
      # natural_max
      #*************
      gene = row[1]
      group = "natural_max"
      position = row[0]
      # 'position', 'gene', 'q_H', 'q_E', 'q_D', 'q_R', 'q_K', 'q_S', 'q_T', 'q_N', 'q_Q', 
      # 'q_A', 'q_V', 'q_L', 'q_I', 'q_M', 'q_F', 'q_Y', 'q_W', 'q_P', 'q_G', 'q_C', 'entropy', 'n_eff'
      res = getMax(row[2:22]) # aa freq
      aa = res[0]
      freq = res[1]

      res_2 = getMax_class(row[24:30]) # aa class freq
      aa_class = res_2[0]
      class_freq = res_2[1]

      writer.writerow([gene, group, position, aa, freq, aa_class, class_freq])

      #*************
      # natural_wt
      #*************
      group = "natural_wt"
      aaList = ['H', 'E', 'D', 'R', 'K', 'S', 'T', 'N', 'Q', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'P', 'G', 'C']
      #print([position, gene, wt[gene]])
      
      try:
        if position in wt[gene]:
          aa = wt[gene][position]
          freq = row[2 + aaList.index(aa)]

          wt_aa_class, wt_class_freq = get_aa_class(aa, class_freqs = row[24:30])
          writer.writerow([gene, group, position, aa, freq, wt_aa_class, wt_class_freq])
        else:
          print(str(gene) + " " + str(position))
      except KeyError:
        print(gene)
      # group = "natural_neff"









