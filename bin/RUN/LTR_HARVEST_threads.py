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
    fastabase=fname.split('/')[-1].replace('.fna','')
    #Function is bieng executed
    os.system(f"gt ltrharvest -index {fname} -maxdistltr {maxdisltr} -mindistltr {mindisltr} -maxlenltr {maxlenltr} -minlenltr {minlenltr} -similar {similar} -seqids yes -gff3 {fname}.harvest.gff3>{fname}.harvest && rm {path}/{fastabase}.fna.des {path}/{fastabase}.fna.esq {path}/{fastabase}.fna.lcp {path}/{fastabase}.fna.llv {path}/{fastabase}.fna.md5 {path}/{fastabase}.fna.prj {path}/{fastabase}.fna.sds {path}/{fastabase}.fna.suf && find . -name '{fastabase}.fna.harvest' -type f -size -570c -delete && find . -name '{fastabase}.fna.harvest.gff3' -type f -size -1c -delete && sed -i '1,12d' {fastabase}.fna.harvest")    
    # set a number of processes to use ncore each
with Pool(int(nprocess)) as p:
    results  = p.map(f, data)
