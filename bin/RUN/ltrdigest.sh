#!/bin/bash
process_id="$1"
userpath="$2"
trna="$3"
LAI="$4"
ltrdigest="$5"
############# use gt to preper the data
cd $userpath/LTRdigest
gt suffixerator -db $userpath/$process_id.fna -indexname $userpath/$process_id.fna -tis -suf -lcp -des -ssp -sds -dna
cp $LAI/$process_id.fna.pass.list.gff3 $ltrdigest/$process_id.fna.pass.list.gff3
sed -E -i 's/\;Classification=\S+//g' $ltrdigest/$process_id.fna.pass.list.gff3
grep -i "##sequence-region" $ltrdigest/$process_id.fna.pass.list.gff3 >$ltrdigest/sequence-region.ids
sed -E -i 's/\s+/\t/g' $ltrdigest/sequence-region.ids
awk -F "\t" '{print $2}' $ltrdigest/sequence-region.ids >$ltrdigest/sequence-region.ids.2
awk '{print $0 "\t" NR-1}' $ltrdigest/sequence-region.ids.2 >$ltrdigest/ids.uniq_with_numbers
sed -i '1,5d' $ltrdigest/$process_id.fna.pass.list.gff3
sed -i 's/Name=/seq_number=/g' $ltrdigest/$process_id.fna.pass.list.gff3
 while read a b
             do
                sed -i "s/$a/seq$b/g" $ltrdigest/$process_id.fna.pass.list.gff3
             done < $ltrdigest/ids.uniq_with_numbers
sed -i 's/seq_number=seq/seq_number=/g' $ltrdigest/$process_id.fna.pass.list.gff3
echo "##gff-version 3" >$ltrdigest/ids.uniq_with_numbers.2
awk -F "\t" '{print "#"$1}' $ltrdigest/ids.uniq_with_numbers >>$ltrdigest/ids.uniq_with_numbers.2
cat ids.uniq_with_numbers.2 $ltrdigest/$process_id.fna.pass.list.gff3 >$ltrdigest/$process_id.fna.pass.list.gff4
gt gff3 -sort $ltrdigest/$process_id.fna.pass.list.gff4 > $ltrdigest/$process_id.fna.pass.list.gff3.sort
sed -i 's/Copia_LTR_retrotransposon/LTR_retrotransposon/g' $ltrdigest/$process_id.fna.pass.list.gff3.sort
sed -i 's/Gypsy_LTR_retrotransposon/LTR_retrotransposon/g' $ltrdigest/$process_id.fna.pass.list.gff3.sort

############## use pass.list.gff3.sort file from LTR_retriever output
gt -j 3 ltrdigest -trnas $userpath/$trna -outfileprefix $process_id  $ltrdigest/$process_id.fna.pass.list.gff3.sort $userpath/$process_id.fna >$ltrdigest/ltrdigest.gff3
