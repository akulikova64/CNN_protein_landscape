import os
import shutil
import csv
import sys

def findMax(ar):
  ind = 2
  m = 0
  for i in range(2,len(ar)):
    if(float(ar[i]) > m):
      ind = i
      m = float(ar[i])
  return [ind, m]
    

#list = os.listdir(path="./aligned_sequences")
#print(list)
#with open('privetik.txt', 'w') as out1:
#shutil.copyfile('privetik.txt','ihha.txt')
#lets do magic epta
'''
path = './data_from_2018_paper/evol_sim_vs_rosetta-master/sequences/designed_sequences_fasta/'

fileList = os.listdir(path)
for file in fileList:
  if file.find('designed') >= 0:
    shutil.copyfile(path + file, './_designed/' + file)
  elif file.find('evolved') >= 0:
    shutil.copyfile(path + file, './_evolved/' + file)
  else:
    print('error in:' + file)
'''
#lets do MORE magic epta
path = './CNN_output_new/'
fileList = os.listdir(path)
aaList = ['H', 'E', 'D',  'R', 'K', 'S', 'T', 'N', 'Q', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'P', 'G', 'C']

with open("cnn_wt_max_freq.csv", "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file)
  writer.writerow(['gene', 'group', 'position', 'aa', 'freq'])

  for file in fileList:
    with open(path + file, 'r') as openedFile:
      lines = [line.rstrip('\n') for line in openedFile]
    #remove header 
    for i in range(len(lines)):
      lines[i] = lines[i].split()
    header = lines[0] 
    del lines[0]
    #pos wt_aa HIS GLU ASP ARG LYS SER THR ASN GLN ALA VAL LEU ILE MET PHE TYR TRP PRO GLY CYS
    # 0     1    2  3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21
    position = 1
    for line in lines:
      # add predicted row
      gene = str(file[0:4]).lower()
      group = "predicted"
      #position = str(line[0])
      res = findMax(line)
      aa = aaList[res[0] - 2]
      freq = str(res[1])
      writer.writerow([gene, group, position, aa, freq])
      # add wt row
      group = "wt"
      aa = aaList[ header.index( line[1] ) - 2 ]
      freq = line[ header.index( line[1] ) ]
      writer.writerow([gene, group, position, aa, freq])

      position += 1
      #sys.exit()
    #print(lines)
    #sys.exit()

