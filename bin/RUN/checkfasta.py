import os
import sys 
import subprocess
from glob import glob
path = sys.argv[1]
fna=False
def is_valid_fasta(dna_sequence):
	for nucleotide in dna_sequence: # we iterate through the chars of our DNA using a foreach
		if nucleotide not in "actgnyACTGNY":
			return False
	return True # if no nucleotide made us return false, it's valid!
 # if no nucleotide made us return false, it's valid!
		
files= glob(f"{path}/*.f*")
if len(files)<3:
    for file in files:
        if '.fna' in file[-5:] or '.fa' in file[-4:] or  '.fasta' in file[-8:]:
            lines=""
            with open(file, "r") as a_file:
                for line in a_file:
                    stripped_line = line.strip()
                    if ">" in stripped_line:
                        continue
                    if len(lines)>=1000:
                        continue
                    lines+=stripped_line
                fna=is_valid_fasta(stripped_line)
else:
    fna = False
   
print(fna)