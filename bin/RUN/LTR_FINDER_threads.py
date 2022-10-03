import os
import glob
import sys

path=sys.argv[1]
nprocess=4
nprocess=sys.argv[2]
LTRFINDER=sys.argv[3]
maxdisltr=sys.argv[4]
mindisltr=sys.argv[5]
maxlenltr=sys.argv[6]
minlenltr=sys.argv[7]
matchpairs=sys.argv[8]
similarFinder=sys.argv[9]
userpath=sys.argv[10]
trna=sys.argv[11]
LAI=sys.argv[12]
FASTA=sys.argv[13]

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

data=glob.glob(f"{path}/*.fna")

from multiprocessing import Pool
def f(fname):
     
    #Function is bieng executed
    os.system(f"{LTRFINDER} {fname} -w 2 -C -D {maxdisltr} -d {mindisltr} -L {maxlenltr} -l {minlenltr} -p {matchpairs} -M {similarFinder} -s {userpath}/{trna} >{fname}.finder")
    #os.system(f"bash threads.sh {x} LTRFINDER")
    

# set a number of processes to use ncore each ### with work as for
with Pool(int(nprocess)) as p:
    results  = p.map(f, data)
