#!/bin/bash
eval "$(conda shell.bash hook)"

TEsorter="$1"
fastafile="$2"
TEsorterdb="$3"
TEsortercov="$4"
TEsortereval="$5"
TEsorterrule="$6"
TEsorterhmm="$7"
RUN="$8"
process_id="$9"
Collected_Files="${10}"
userpath="${11}"
threads="${12}"
mkdir -p $TEsorter/tmp
tmp=$TEsorter/tmp
# echo "$userpath"
conda activate MegaLTR
TEsorter $fastafile -db $TEsorterdb -cov $TEsortercov -eval $TEsortereval  -rule $TEsorterrule --hmm-database $TEsorterhmm -tmp $tmp -p $threads 2>/dev/null