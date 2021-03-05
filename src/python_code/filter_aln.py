# this program sorts alignments 
import os
import sys
from shutil import copyfile

perc_sims = ["20","40","60","80", "100"]
min_lines = 50 * 2

for perc_sim in perc_sims:
  input_path = "../../data/PSICOV/aln_" + perc_sim + "/"
  output_path = "../../data/PSICOV/aln_filtered/aln_" + perc_sim + "/"

  protein_list = os.listdir(input_path)

  for protein in protein_list:
    with open(input_path + protein, 'r') as file:
      lines = [line.rstrip('\n') for line in file]
    if len(lines) >= min_lines:
      copyfile(input_path + protein, output_path + protein)
  #with open(new_path, "w") as file:
