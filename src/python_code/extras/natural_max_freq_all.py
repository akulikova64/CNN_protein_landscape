import os
import shutil
import csv
import sys

def getMax(list, gene, position):
  aaList = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]  
  ind = 0
  max = 0

  for i in range(0,len(list)):
    if(float(list[i]) > max):
      ind = i
      max = float(list[i])

  # now just looking for ties with the detected maximum
  aa_tie_list = []
  tie_count = 0

  for i in range(0, len(list)):
    current = float(list[i])
    if(round(current, 5) == round(max, 5)):
      tie_count += 1
      aa_tie_list.append(str(aaList[i]))

  '''if(tie_count > 1):

    print(str(gene), "position:", str(position))
    print("max_freq:", str(round(max, 5)))
    print("tie count:", str(tie_count - 1)) # not counting the first highest amino acid freq found
    print("tied aa's:", str(aa_tie_list))
    print()'''
  
  return [aaList[ind], max]

def getMax_class(list):
  class_list = ['unique', 'aliphatic', 'small_polar', 'negative', 'positive', 'aromatic']
  ind = 0
  max = 0
  for i in range(0, len(list)):
    if(float(list[i]) > max):
      ind = i
      max = float(list[i])

  return [class_list[ind], max]

def get_aa_class(wt_aa, class_freqs):
  class_list = ['unique', 'aliphatic', 'small_polar', 'negative', 'positive', 'aromatic']

  unique = ["P", "G"]
  aliphatic = ["M", "L", "I", "V", "A"]
  small_polar = ["C", "S", "T", "N", "Q"]
  negative = ["D", "E"]
  positive = ["R", "K"]
  aromatic = ["H", "Y", "F", "W"]

  # get aa_class of the wt (from "wt_aa" parameter)  
  if wt_aa in unique:
    aa_class = "unique"
  if wt_aa in aliphatic:
    aa_class = "aliphatic" 
  if wt_aa in small_polar:
    aa_class = "small_polar"
  if wt_aa in negative:
    aa_class = "negative"
  if wt_aa in positive:
    aa_class = "positive"
  if wt_aa in aromatic:
    aa_class = "aromatic"

  ind = class_list.index(aa_class)
  class_freq = class_freqs[ind]

  return aa_class, class_freq

#--------- main ------------------------------------------
box_size_list = ["12", "20", "30", "40"]

for box_size in box_size_list:

  input_path_1 = "../../../data/PSICOV_box_" + box_size + "/output/cnn_wt_max_freq.csv"
  input_path_2 = "../../../output/output_PSICOV/stats_align_all.csv"
  output_path = "../../../data/PSICOV_box_" + box_size + "/output/natural_max_freq_files/natural_max_freq_all.csv" 

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
        
        # 'position,gene,q_A,q_R,q_N,q_D,q_C,q_Q,q_E,q_G,q_H,q_I,q_L,q_K,q_M,q_F,q_P,q_S,q_T,q_W,q_Y,q_V,entropy,n_eff,q_aliphatic,q_polar,q_positive,q_negative,q_aromatic,q_proline,entropy_class,n_eff_class

        res = getMax(row[2:22], gene, position) # aa freq
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
        aaList = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
        #print([position, gene, wt[gene]])
        
        try:
          if position in wt[gene]:
            aa = wt[gene][position]
            freq = row[2 + aaList.index(aa)]

            wt_aa_class, wt_class_freq = get_aa_class(aa, class_freqs = row[24:30])
            writer.writerow([gene, group, position, aa, freq, wt_aa_class, wt_class_freq])
          else:
            continue
            #print(str(gene), str(position), "position not in both cnn and MSA")
        except KeyError:
          continue
          #print(gene, "gene not found")
        # group = "natural_neff"


#supplementary info: counting the number ot ties:

def get_ties(list, gene, position):
  aaList = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]  
  ind = 0
  max = 0

  for i in range(0,len(list)):
    if(float(list[i]) > max):
      ind = i
      max = float(list[i])

  # now just looking for ties with the detected maximum
  aa_tie_list = []
  tie_count = 0

  for i in range(0, len(list)):
    current = float(list[i])
    if(round(current, 5) == round(max, 5)):
      tie_count += 1
      aa_tie_list.append(str(aaList[i]))
  
  return [gene, position, aa_tie_list, (tie_count - 1), round(max, 5)]

for box_size in box_size_list:
  output_txt = "../../../data/PSICOV_box_" + box_size + "/output/consensus_ties.csv" 
  with open(output_txt, "w", newline='\n', encoding='utf-8') as CSV_file:
    writer = csv.writer(CSV_file)
    writer.writerow(['gene', 'position', 'aa_list', 'tie_count', 'freq'])

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
        position = row[0]
        
        res = get_ties(row[2:22], gene, position) # aa freq
        writer.writerow([res[0], res[1], res[2], res[3], res[4]])