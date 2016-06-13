# Project 2 proteomics : dealing with the fasta files:
import pandas as pd
import numpy as np
from sys import argv

amino_acid_masses = { 'A': 71.0788,'R': 156.1875,'N': 114.1038,'D': 115.0886,'C': 103.1388,'E': 129.1155,
    'Q': 128.1307,'G': 57.0519,'H': 137.1411,'I': 113.1594,'L': 113.1594,'K': 128.1741,'M': 131.1926,'F': 147.1766,
    'P': 97.1167,'S': 87.0782,'T': 101.1051,'W': 186.2132,'Y': 163.1760,'V': 99.1326,
    'X': 136.9008-18,'Z': 146.6375-18,'B': 132.6108-18 }

# Opening and reading the fasta file
file_open = open(argv[1], 'r').read()
file_separe = file_open.split('>') 
header = []
file_separe.remove('')
sequence_dict = {}
for entry in file_separe:
		seq = entry.splitlines()
		if seq[0] not in sequence_dict:
			sequence_dict[seq[0]] = (seq[1])

# Spliting the proteins imitating the trypsine enzyme
def trypsin_cut(protein):
	cleaved = []
	prev_i = 0
	for i in range(len(protein)-1):
		if ((protein[i]=='K' or protein[i]=='R') and (protein[i+1]!='P')):
			cleaved.append(protein[prev_i:i+1])
			prev_i = i+1
	
	return cleaved

def compute_mass(s):
	total = 0
	for aa in s:
		total += amino_acid_masses[aa]
	return total

# Calculating the fragments of the entire data:
mass = []
mass_temp = 0
protname = []
protmass = []
which = []
seq = []
for key, values in sequence_dict.items():
	fragments = trypsin_cut(values)
	for j in fragments:
		mass_temp = compute_mass(j)
		if mass_temp > 700 and mass_temp < 3000:
			mass.append(compute_mass(j))
			protmass.append(compute_mass(values))
			protname.append(key)
			seq.append(j)

# Saving results in a Pandas dataframe format
result = pd.DataFrame({'Seqs': seq, 'Protein': protname, 'Peptide Mass': mass, 'Protein Mass': protmass})
#description = result.head()
result.to_csv('library.csv') 

print result[:30]