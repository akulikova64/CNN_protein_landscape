import csv
import os
# get the average % homology of all alignments

def get_percent(aln_seq, ref_seq):
  count_similar = 0
  for i, j in zip(aln_seq, ref_seq):
    assert j != "-"
    if i == j:
      count_similar += 1s

  similarity = float(count_similar/len(ref_seq))
  return similarity

#------------------------------------main-------------------------------------------------
input_path = "../../data/PSICOV/aln_fasta/"
output_path = "../../data/PSICOV/" 

protein_list = os.listdir(input_path)

for protein in protein_list:
  with open(output_path + protein, 'w') as csv_file:
    records = list(SeqIO.parse(input_path + protein, "fasta"))
    ref_seq = records[0].seq
    #file.write(">reference\n") # we do not want to include the reference in the calculations
    #file.write(str(ref_seq) + "\n")

    for i in range(1, len(records)):
      aln_seq = records[i].seq
      keep = compare(str(aln_seq), str(ref_seq))

      if keep:
        file.write(">" + str(i) + "\n")
        file.write(str(aln_seq) + "\n")

      
  