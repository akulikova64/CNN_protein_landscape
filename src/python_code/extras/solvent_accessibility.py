# this program gets solven accessibility values.

import csv
import math
import os
import sys

# extracting the solvent accessibility from the CIF output files. 
aa_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR','TRP', 'TYR', 'VAL']
aaCodes = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T','TRP':'W', 'TYR':'Y', 'VAL':'V'}

input_path = "../../../data/PSICOV/psicov_CIF/" #input multiple sequence aligments path
output_path = "../../../output/output_PSICOV/solvent_accessibility.csv"

#list of CIF files:
protein_list = os.listdir(input_path)

with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:
  writer = csv.writer(CSV_file) 

  writer.writerow(['gene', 'position', 'aa', 'SASA_abs_total', 'SASA_rel_total'])

  for protein in protein_list:
    with open(input_path + protein, 'r') as file:
    
      gene = protein[0:4]
      detected = False
      for line in file:
        line = line.strip('\n')
        line = line.split()

        if line and line[0] == "_freeSASA_rsa.rel_polar":
          detected = True
          continue
          
        if detected:
          residue = line[2]

          if residue in aa_list:
            position = line[1]
            aa = aaCodes[residue]
            sasa_abs = line[3]
            sasa_rel = line[4]

            writer.writerow([str(gene), str(position), str(aa), str(sasa_abs), str(sasa_rel)]) 

print("Saved CSV to " + output_path)