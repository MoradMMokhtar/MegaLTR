#!/bin/bash
process_id="$1"
userpath="$2"
trna="$3"
LAI="$4"
ltrdigest="$5"
############# use gt to preper the data
cd $userpath/LTRdigest
gt suffixerator -db $userpath/$process_id.fna -indexname $userpath/$process_id.fna -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index $userpath/$process_id.fna -gff3 $LAI/$process_id.fna.pass.list.gff3 >$ltrdigest/for_dlete2
gt gff3 -sort $LAI/$process_id.fna.pass.list.gff3 > $ltrdigest/$process_id.fna.pass.list.gff3.sort
############## use pass.list.gff3 file from LTR_retriever output
gt -j 3 ltrdigest -trnas $userpath/$trna -outfileprefix $process_id  $ltrdigest/$process_id.fna.pass.list.gff3.sort $userpath/$process_id.fna >$ltrdigest/ltrdigest.gff3
