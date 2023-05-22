#!/bin/bash
process_id="$1"
userpath="$2"
trna="$3"
LAI="$4"
ltrdigest="$5"
############# use gt to preper the data
   cd $userpath/LTRdigest
   ln -s $userpath/$process_id.fna $userpath/LTRdigest/$process_id.fna > /dev/null 2>/dev/null
   gt suffixerator -db $process_id.fna -indexname $process_id.fna -tis -suf -lcp -des -ssp -sds -dna
   
   if [ -f "$LAI/$process_id.fna.pass.list.gff3" ]; then 
   cp $LAI/$process_id.fna.pass.list.gff3 $ltrdigest ##GFF3 format for intact LTR-RTs
   fi
   if [ -f "$LAI/$process_id.fna.mod.pass.list.gff3" ]; then 
   cp $LAI/$process_id.fna.mod.pass.list.gff3 $ltrdigest/$process_id.fna.pass.list.gff3 ##GFF3 format for intact LTR-RTs
   fi
   
   sed -E -i 's/\;Classification=\S+//g' $ltrdigest/$process_id.fna.pass.list.gff3
   grep -i "##sequence-region" $ltrdigest/$process_id.fna.pass.list.gff3 >$ltrdigest/sequence-region.ids
   sed -E -i 's/\s+/\t/g' $ltrdigest/sequence-region.ids
   awk -F "\t" '{print $2}' $ltrdigest/sequence-region.ids >$ltrdigest/sequence-region.ids.2
   cp $ltrdigest/$process_id.fna.des $ltrdigest/$process_id.fna.ids
   sed -i '$d' $ltrdigest/$process_id.fna.ids
   sed -i '/^$/d' $ltrdigest/$process_id.fna.ids
   awk '{print $0 "\t" NR-1}' $ltrdigest/$process_id.fna.ids >$ltrdigest/ids.uniq_with_numbers
   sed -i '1,5d' $ltrdigest/$process_id.fna.pass.list.gff3
   sed -i 's/Name=/seq_number=/g' $ltrdigest/$process_id.fna.pass.list.gff3
    while read a b
                do
                   sed -i "s/$a/seq$b/g" $ltrdigest/$process_id.fna.pass.list.gff3
                done < $ltrdigest/ids.uniq_with_numbers
   sed -i 's/seq_number=seq/seq_number=/g' $ltrdigest/$process_id.fna.pass.list.gff3
   echo "##gff-version 3" >$ltrdigest/ids.uniq_with_numbers.2
   awk -F "\t" '{print "#"$1}' $ltrdigest/sequence-region.ids.2 >>$ltrdigest/ids.uniq_with_numbers.2
   cat $ltrdigest/ids.uniq_with_numbers.2 $ltrdigest/$process_id.fna.pass.list.gff3 >$ltrdigest/$process_id.fna.pass.list.gff4
   gt gff3 -sort $ltrdigest/$process_id.fna.pass.list.gff4 > $ltrdigest/$process_id.fna.pass.list.gff3.sort
   sed -i 's/Copia_LTR_retrotransposon/LTR_retrotransposon/g' $ltrdigest/$process_id.fna.pass.list.gff3.sort
   sed -i 's/Gypsy_LTR_retrotransposon/LTR_retrotransposon/g' $ltrdigest/$process_id.fna.pass.list.gff3.sort

############## use pass.list.gff3.sort file from LTR_retriever output
   gt -j 3 ltrdigest -trnas $userpath/$trna -outfileprefix $process_id  $ltrdigest/$process_id.fna.pass.list.gff3.sort $ltrdigest/$process_id.fna >$ltrdigest/ltrdigest.gff3
