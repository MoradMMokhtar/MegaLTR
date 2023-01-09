import os
import glob
import sys

# print ("this parameters 1", sys.argv)
fastafile=sys.argv[1]
LTRFINDEr=sys.argv[2]  
maxdisltr=sys.argv[3]
mindisltr=sys.argv[4]
maxlenltr=sys.argv[5]
minlenltr=sys.argv[6]
matchpairs=sys.argv[7]
similarFinder=sys.argv[8]
userpath=sys.argv[9]
trna=sys.argv[10]

script=sys.argv[11]
fasta=sys.argv[12]
LTRHARVESt=sys.argv[13] 
LAi=sys.argv[14]
similar=sys.argv[15]
procces_id=sys.argv[16]
nprocess=sys.argv[17]
# set number of CPUs to run on
ncore = f"{nprocess}"

# ncore = "4"

# set env variables
# have to set these before importing numpy
os.environ["OMP_NUM_THREADS"] = ncore
os.environ["OPENBLAS_NUM_THREADS"] = ncore
os.environ["MKL_NUM_THREADS"] = ncore
os.environ["VECLIB_MAXIMUM_THREADS"] = ncore
os.environ["NUMEXPR_NUM_THREADS"] = ncore

# data=glob.glob(f"{path}/*.fna")


# print(ids)

#read the fasta file 
sequences = []
ids = []
descr = None
i = 0
with open(fastafile) as file:
	line = file.readline()[:-1] # always trim newline
	while line:
		if line[0] == '>':
			if descr: # any sequence found yet?
				sequences.append((descr, seq,i))
				i+=1
			descr = str(line[1:].split('>')).replace("]","").replace("[","").replace("'","")
			seq = '' # start a new sequence
		else:
			seq += line
		line = file.readline()[:-1]
	sequences.append((descr, seq,i))
	i+=1


from multiprocessing import Pool
def f(val):
	chr_id=val[2]
	chr_file=f"{fasta}/{val[0].split()[0]}.fna"
	with open(chr_file,"w") as file:
		header_fasta=f"{val[0]}"
		chr_fasta=val[1]
		file.write(f">{header_fasta}\n")
		file.write(f"{chr_fasta}\n")
	

    #Function is bieng executed
    #os.system(f"{LTRFINDER} {fname} -w 2 -C -D {maxdisltr} -d {mindisltr} -L {maxlenltr} -l {minlenltr} -p {matchpairs} -M {similarFinder} -s {userpath}/{trna} >{fname}.finder && find . -name '{fastabase}.fna.finder' -type f -size -398c -delete && cat {fname}.finder >>{FASTA}/all.fna.finder")
	os.system(f"bash {script}/LTR_paralelle.sh  {chr_file} {LTRFINDEr} {maxdisltr} {mindisltr} {maxlenltr} {minlenltr} {matchpairs} {similarFinder} {similar} {userpath} {trna} {fasta} {LTRHARVESt} {LAi} {procces_id} {script} {chr_id} {val[0].split()[0]}")


    #os.system(f"bash threads.sh {x} LTRFINDER")
    
# set a number of processes to use ncore each ### with work as for
with Pool(int(nprocess)) as p:
    results  = p.map(f, sequences)
