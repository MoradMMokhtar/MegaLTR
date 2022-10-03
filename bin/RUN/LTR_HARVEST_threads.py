import os
import glob
import sys

path=sys.argv[1]
nprocess=4
nprocess=sys.argv[2]
maxdisltr=sys.argv[3]
mindisltr=sys.argv[4]
maxlenltr=sys.argv[5]
minlenltr=sys.argv[6]
matchpairs=sys.argv[7]
similar=sys.argv[8]
LTRHARVEST=sys.argv[9]
# set number of CPUs to run on
ncore = f"{nprocess}"
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
    os.system(f"gt ltrharvest -index {fname} -maxdistltr {maxdisltr} -mindistltr {mindisltr} -maxlenltr {maxlenltr} -minlenltr {minlenltr} -similar {similar} -seqids yes -gff3 {fname}.harvest.gff3>{fname}.harvest")    
    # set a number of processes to use ncore each
with Pool(int(nprocess)) as p:
    results  = p.map(f, data)
