import json
import sys

file_num = ["", "2", "3"]
input_path = "../../../../training_data/snapshots_cif_uniq_modCIF_chain_dataset_no_psicov_50_seq_sim"
output_path = "../../../../training_data/pdbs.txt"

pdbs = []

for num in file_num:
  with open(input_path + num + ".txt", "r") as json_file:
    for line in json_file:
      items = json.loads(line)
      pdb = items["pdb_code"]
      if pdb not in pdbs:
        pdbs.append(pdb)
    print("got all pdbs, making text file")

with open(output_path, "w") as pdb_file:
  for pdb in pdbs:
    pdb_file.write(pdb + "\n")
    