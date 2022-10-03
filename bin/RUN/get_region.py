import sys
fasta_path=sys.argv[1]
name=sys.argv[2]
output=sys.argv[3]
data=[]
#print(name,fasta_path,output)
#parametre example  python3 get_region.py files/Arabidopsis_thaliana.fna LTR_Table_TE_TEsorter_EDTA_Digest outputfile
with open(f"{name}",'r') as f:
    ftmp =f.readlines()
    for fm in ftmp:
        file=fm.split()
        data.append((file[0],file[1],file[5],file[6],file[8],file[9]))#,file[3]))
sequences = {}
descr = None
with open(f"{fasta_path}","r") as file:
    line = file.readline()[:-1] # always trim newline
    while line:
        if line[0] == '>':
            if descr: # any sequence found yet?
                 sequences[descr.split()[0][1:]]= seq
            descr = str(line)

            seq = '' # start a new sequence
        else: 
            seq += line
        line = file.readline()[:-1]
    sequences[descr.split()[0][1:]]= seq
for line in data:
    with open(f"{output}/{line[0]}.fasta",'w') as file:
        if line[2].isnumeric() and line[3].isnumeric() and line[5].isnumeric() and line[4].isnumeric():
            file.write(">5'\n")
            file.write(str(sequences[line[1]][int(line[2]):int(line[3])])+"\n")
           
            file.write(">3'\n")
            file.write(str(sequences[line[1]][int(line[4]):int(line[5])])+"\n")
           
