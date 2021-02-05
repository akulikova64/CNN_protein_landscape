from Bio import SeqIO

# This code samples the sequences PSICOV alignments by % similarity to the sequences in the PSICOV pdb structures.

#Plan:
#load the seqs from the CNN output (wt aa) into a dictionary with the protein name. 
#open and parse each of  the fasta files 
# make sure that name of fasta file matches with name the pdb sequences. 

#loading a fasta aln file:

#THRESHOLD = 0.20
THRESHOLD = 1.00

def compare(aln_seq, ref_seq):
  count_similar = 0
  for i, j in zip(aln_seq, ref_seq):
    if i == j:
      count_similar += 1
    else:
      print("aln:", i, "ref", j)
  
  print(count_similar)
  print(len(aln_seq))
  max_limit = THRESHOLD + 0.05
  min_limit = THRESHOLD - 0.05

  similarity = count_similar/len(aln_seq)

  if similarity >= min_limit and similarity <= max_limit:
    #print(similarity)
    return True
  else:
    #print(similarity)
    return False

def get_reference_from_CNN_data():
  aaCodes = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T','TRP':'W', 'TYR':'Y', 'VAL':'V'}
  with open("../data/PSICOV/PSICOV_CNN_output/1a3a_final_tot.csv", 'r') as file:
    ref_seq = ""
    file.readline()
    for line in file:
      line = line.rstrip('\n')
      line = line.rsplit(",")
      wt_aa = line[5]
      ref_seq += aaCodes[wt_aa]

  return ref_seq

CNN_reference = get_reference_from_CNN_data()
records = list(SeqIO.parse("../data/PSICOV/aln_fasta/1a3aA.fasta", "fasta"))
#ref_seq = get_reference()
ref_seq = records[0].seq

#print(CNN_reference[0:40])
#print(ref_seq[0:40])
print(compare(CNN_reference, ref_seq))

'''
with open("../data/PSICOV/aln_20/1a3a_aln_20.fasta", 'w') as file:
  

  records = list(SeqIO.parse("../data/PSICOV/aln_fasta/1a3aA.fasta", "fasta"))
  ref_seq = get_reference()

  for i in range(len(records)):
    aln_seq = records[i].seq
    keep = compare(aln_seq, ref_seq)

    if keep:
      # then write the aln_seq to a new file and folder (aln_20)'''
      

