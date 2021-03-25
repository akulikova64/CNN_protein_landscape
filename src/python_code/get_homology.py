from Bio import SeqIO
import csv
import os
# get the average % homology of all alignments

def get_percent(aln_seq, ref_seq):
  count_similar = 0
  for i, j in zip(aln_seq, ref_seq):
    assert j != "-"
    if i == j:
      count_similar += 1

  similarity = float(count_similar/len(ref_seq))
  return similarity

#---------------main-------------------------
input_path = "../../data/PSICOV/aln_fasta/"
output_path = "../../data/PSICOV/" 


with open(output_path + "percent_homol.csv", 'w', newline='\n', encoding='utf-8') as csv_file:
  writer = csv.writer(csv_file)
  writer.writerow(['gene', 'percent_homol'])

  protein_list = os.listdir(input_path)

  for protein in protein_list:
    records = list(SeqIO.parse(input_path + protein, "fasta"))
    ref_seq = records[0].seq

    percent_sum = 0
    for i in range(1, len(records)):
      aln_seq = records[i].seq
      percent_sum += get_percent(str(aln_seq), str(ref_seq))

    percent = percent_sum/len(records)
    writer.writerow([protein[0:4], percent])

      
      
  