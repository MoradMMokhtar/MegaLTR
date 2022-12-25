import os
import glob
import sys

FASTAfile=sys.argv[1]
LTRfiles=sys.argv[2]
outdir=sys.argv[3]
nprocess=sys.argv[4]
extractseq=sys.argv[5]
# nprocess=6
# set number of CPUs to run on
ncore = f"{nprocess}"
# ncore = "4"
os.environ["OMP_NUM_THREADS"] = ncore
os.environ["OPENBLAS_NUM_THREADS"] = ncore
os.environ["MKL_NUM_THREADS"] = ncore
os.environ["VECLIB_MAXIMUM_THREADS"] = ncore
os.environ["NUMEXPR_NUM_THREADS"] = ncore

data=glob.glob(f"{LTRfiles}/*")

from multiprocessing import Pool
def f(fname):
     
    #Function is bieng executed
    os.system(f"perl {extractseq} {FASTAfile} {fname} {outdir}/LTR-RT_Sequence.fa")

# set a number of processes to use ncore each ### with work as for
with Pool(int(nprocess)) as p:
    results  = p.map(f, data)
