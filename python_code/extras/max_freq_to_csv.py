import os
import shutil
import csv
import sys

def findMax(line):
  ind = 2
  m = 0
  for i in range(2, len(line)):
    if(float(line[i]) > m):
      ind = i
      m = float(line[i])
  return [ind, m]

def findMax_class(class_dict):
  aa_class_name = ""
  max_class_freq = 0 # max
  for aa_class in class_dict:
    if class_dict[aa_class] > max_class_freq:
      aa_class_name = aa_class
      max_class_freq = class_dict[aa_class]

  return aa_class_name, max_class_freq
    
def get_class_freq(line):

  # old: aaList = ['H', 'E', 'D',  'R', 'K', 'S', 'T', 'N', 'Q', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'P', 'G', 'C']
  aaList = ['A', 'R', 'N',  'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

  aliphatic = ["G", "A", "V", "L", "M", "I"]
  polar = ["S", "T", "C", "N", "Q"]
  positive = ["K", "R", "H"]
  negative = ["D", "E"]
  aromatic = ["F", "Y", "W"]
  proline = "P"

  class_dict = {"aliphatic":0, "polar":0, "positive":0, "negative":0, "aromatic":0, "proline":0}
  for i, aa in enumerate(aaList):
    if aa in aliphatic:
      class_dict["aliphatic"] += float(line[i+10])
    if aa in polar:
      class_dict["polar"] += float(line[i+10]) 
    if aa in positive:
      class_dict["positive"] += float(line[i+10])
    if aa in negative:
      class_dict["negative"] += float(line[i+10])
    if aa in aromatic:
      class_dict["aromatic"] += float(line[i+10])
    if aa == proline:
      class_dict["proline"] += float(line[i+10])

  aa_class, class_freq = findMax_class(class_dict)

  return aa_class, class_freq

def get_class_freq_2(wt_aa, line):
  # old: aaList = ['H', 'E', 'D',  'R', 'K', 'S', 'T', 'N', 'Q', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'P', 'G', 'C']
  aaList = ['A', 'R', 'N',  'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

  aliphatic = ["G", "A", "V", "L", "M", "I"]
  polar = ["S", "T", "C", "N", "Q"]
  positive = ["K", "R", "H"]
  negative = ["D", "E"]
  aromatic = ["F", "Y", "W"]
  proline = "P"

  class_dict = {"aliphatic":0, "polar":0, "positive":0, "negative":0, "aromatic":0, "proline":0}
  for i, aa in enumerate(aaList):
    if aa in aliphatic:
      class_dict["aliphatic"] += float(line[i+10])
    if aa in polar:
      class_dict["polar"] += float(line[i+10]) 
    if aa in positive:
      class_dict["positive"] += float(line[i+10])
    if aa in negative:
      class_dict["negative"] += float(line[i+10])
    if aa in aromatic:
      class_dict["aromatic"] += float(line[i+10])
    if aa == proline:
      class_dict["proline"] += float(line[i+10])

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

  class_freq = class_dict[aa_class]

  return aa_class, class_freq

#--------------- main ----------------------------
input_path = '../../data/PSICOV/PSICOV_CNN_output/'
output_path = "../../data/PSICOV/output/cnn_wt_max_freq.csv"

fileList = os.listdir(input_path)
# old: aaList = ['H', 'E', 'D',  'R', 'K', 'S', 'T', 'N', 'Q', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'P', 'G', 'C']
aaList = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
aaCodes = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T','TRP':'W', 'TYR':'Y', 'VAL':'V'}

with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file)
  writer.writerow(['gene', 'group', 'position', 'aa', 'freq', 'aa_class', 'class_freq'])

  for file in fileList:
    with open(input_path + file, 'r') as openedFile: # input is CNN output file
      lines = [line.rstrip('\n') for line in openedFile]
    #remove header 
    for i in range(len(lines)):
      lines[i] = lines[i].split(",")
    header = lines[0] 
    del lines[0]
    # old: structure of "line" in lines:
    #pos wt_aa HIS GLU ASP ARG LYS SER THR ASN GLN ALA VAL LEU ILE MET PHE TYR TRP PRO GLY CYS
    # 0     1    2  3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21

    # new (correct) structure of "line" in lines:
    #pos, aa_id, pdb_id, chain_id, pos, wtAA, prAA, wt_prob, pred_prob, avg_log_ratio, prALA, prARG, prASN, prASP, prCYS,prGLN,prGLU,prGLY,prHIS,prILE,prLEU,prLYS,prMET,prPHE,prPRO,prSER,prTHR,prTRP,prTYR,prVAL,prHydrophobic,prAromatic,prPolarUncharged,prCationic,prAnionic,prCharged,prSmall,prSulfur,prAcyl,prAlcohol
    # 0     1      2        3       4    5     6      7          8            9         10      11     12    13     14    15    16    17   18    19    20    21     22     23    24   25    26    27    28    29      30             31         32               33        34
    position = 1
    for line in lines:

      gene = str(file[0:4]).lower()

      # add predicted row
      group = "predicted"
      #position = str(line[0])
      #res = findMax(line)
      #aa = aaList[res[0] - 2]
      aa = aaCodes[str(line[6])]
      freq = str(line[8])
      aa_class, class_freq = get_class_freq(line)
      writer.writerow([gene, group, position, aa, freq, aa_class, class_freq])

      # add wt row
      group = "wt"
      # old : wt_aa = aaList[ header.index( line[1] ) - 2 ]
      wt_aa = aaCodes[str(line[5])]
      # old : freq = line[ header.index( line[1] ) ]
      freq = str(line[7])
      wt_aa_class, wt_class_freq = get_class_freq_2(wt_aa, line)
      writer.writerow([gene, group, position, wt_aa, freq,  wt_aa_class, wt_class_freq])

      position += 1
      #sys.exit()
    #print(lines)
    #sys.exit()

