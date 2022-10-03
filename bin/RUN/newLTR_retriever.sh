#!/bin/bash
genome="$1"
userpath="$2"
LAI="$3"
minlen="$4"
process_id="$5"
Collected_Files="$6"
threads="$7"
LTRretriever="$8"
eval "$(conda shell.bash hook)"
cd $LAI
conda activate MegaLTR

#$LTR_retriever -genome $genome -inharvest $LAI/all.harvest.finder.combine -threads $threads -minlen $minlen >>$userpath/screen.txt
$LTRretriever -genome $genome -inharvest $LAI/all.harvest.finder.combine -threads $threads -minlen $minlen >>$userpath/screen.txt
cp $LAI/$process_id.fna.pass.list $Collected_Files ##All LTR-RTs
cp $LAI/$process_id.fna.nmtf.pass.list $Collected_Files ##Non-TGCA LTR-RTs
cp $LAI/$process_id.fna.pass.list.gff3 $Collected_Files ##GFF3 format for intact LTR-RTs
cp $LAI/$process_id.fna.LTRlib.redundant.fa $Collected_Files ##All LTR-RTs with redundancy
cp $LAI/$process_id.fna.LTRlib.fa $Collected_Files ##All non-redundant LTR-RTs
cp $LAI/$process_id.fna.LTR.gff3 $Collected_Files ##GFF3 format
cp $LAI/$process_id.fna.out.fam.size.list $Collected_Files ##LTR family summary
cp $LAI/$process_id.fna.out.superfam.size.list $Collected_Files ##LTR superfamily summary
cp $LAI/$process_id.fna.defalse $Collected_Files
cp $LAI/$process_id.fna.masked $Collected_Files
cp $LAI/$process_id.fna.out $Collected_Files
cp $LAI/$process_id.fna.tbl $Collected_Files
