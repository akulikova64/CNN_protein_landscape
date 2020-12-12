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


#--------- main ------------------------------------------
input_path_1 = "../../data/output/cnn_wt_max_freq.csv"
input_path_2 = "../../data/output/stats_align_all.csv"
output_path = "../../data/output/natural_max_freq.csv" 

wt = {}
with open(input_path_1, "r", newline='\n', encoding='utf-8') as CSV_file:
  csv_reader = csv.reader(CSV_file, delimiter=',')
  # gene, group, position, aa, freq
  #  0     1        2      3    4
  gene = ""
  for row in csv_reader:
    if row[1] == "wt":
      if row[0] != gene:
        wt[row[0]] = {}
        gene = row[0]
      wt[row[0]][row[2]] = row[3]

with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file)
  writer.writerow(['gene', 'group', 'position', 'aa', 'freq'])

  with open(input_path_2, "r", newline='\n', encoding='utf-8') as CSV_file:
    csv_reader = csv.reader(CSV_file, delimiter=',')
    for row in csv_reader:
      if row[0] == 'position':
        continue

      #*************
      # natural_max
      #*************
      gene = row[1]
      group = "natural_max"
      position = row[0]
      # 'position', 'gene', 'q_H', 'q_E', 'q_D', 'q_R', 'q_K', 'q_S', 'q_T', 'q_N', 'q_Q', 
      # 'q_A', 'q_V', 'q_L', 'q_I', 'q_M', 'q_F', 'q_Y', 'q_W', 'q_P', 'q_G', 'q_C', 'entropy', 'n_eff'
      res = getMax(row[2:22])
      aa = res[0]
      freq = res[1]
      writer.writerow([gene, group, position, aa, freq])

      #*************
      # natural_wt
      #*************
      group = "natural_wt"
      aaList = ['H', 'E', 'D', 'R', 'K', 'S', 'T', 'N', 'Q', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'P', 'G', 'C']
      #print([position, gene, wt[gene]])
      if position in wt[gene]:
        aa = wt[gene][position]
        freq = row[2 + aaList.index(aa)]
        writer.writerow([gene, group, position, aa, freq])
      else:
        print(str(gene) + " " + str(position))
      # group = "natural_neff"









