from Bio import SeqIO
import os
import sys


# this code aligns a protein sequence alignment reference to the wild type sequence in the CNN output

def get_reference_from_CNN_data(name):
  aaCodes = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T','TRP':'W', 'TYR':'Y', 'VAL':'V'}
  ref_seq = ""

  try:
    file = open("../../data/PSICOV/PSICOV_CNN_output/" + name + "_final_tot.csv", 'r')
  except FileNotFoundError:
    try:
      file = open("../../data/PSICOV/PSICOV_CNN_output/" + name + ".csv", 'r')
    except FileNotFoundError:
      return "missing"

  file.readline()

  for line in file:
    line = line.rstrip('\n')
    line = line.rsplit(",")
    wt_aa = line[5]
    ref_seq += aaCodes[wt_aa]

  file.close()

  return str(ref_seq)

def compare(aln_ref, CNN_ref):
  length = len(aln_ref)
  CNN_ref = CNN_ref[0:length]
  count_similar = 0

  for i, j in zip(aln_ref, CNN_ref):
    if i == j:
      count_similar += 1

  similarity = float(count_similar/length)

  return similarity

input_path = "../../data/PSICOV/aln_fasta/"
protein_list = os.listdir(input_path)

unaligned_count = 0
aligned_count = 0
missing_count = 0
unaligned_list = []

for protein in protein_list:
  records = list(SeqIO.parse(input_path + protein, "fasta"))
  aln_ref = records[0].seq

  name = protein[0:4]
  CNN_ref = get_reference_from_CNN_data(name)

  if CNN_ref != "missing":
    similarity = compare(aln_ref, CNN_ref)
    if similarity != 1.0:
      unaligned_count += 1
      unaligned_list.append(str(name))
      #print(name, similarity)
      #print(aln_ref)
      #print(CNN_ref)
      #print()
    else:
      aligned_count += 1
  else:
    print(name, "missing")
    missing_count += 1

print()
print("missing:", missing_count)
print("unaligned:", unaligned_count)
print("aligned:", aligned_count)
print()
print("list unaligned:", unaligned_list)