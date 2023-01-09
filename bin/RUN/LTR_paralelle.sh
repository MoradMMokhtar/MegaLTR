#!/bin/bash
fastafile=${1}
LTRFINDER=${2}
maxdisltr=${3}
mindisltr=${4}
maxlenltr=${5}
minlenltr=${6}
matchpairs=${7}
similarFinder=${8}
similar=${9}
out=${10}
trna=${11}
fasta=${12}
LTRHARVEST=${13}
LAI=${14}
procces_id=${15}
RUN=${16}
ID_chromo=${17}
NC_chromo=${18}
eval "$(conda shell.bash hook)"
conda activate MegaLTR

  gt suffixerator -db $fastafile -indexname $fastafile -tis -suf -lcp -des -ssp -sds -dna

  gt ltrharvest -index $fastafile -maxdistltr $maxdisltr -mindistltr $mindisltr -maxlenltr $maxlenltr -minlenltr $minlenltr -similar $similar -seqids yes -gff3 $fastafile.harvest.gff3> $fastafile.harvest
  
  dirname="$(dirname "${fastafile}")"
  fastaname=$(basename $fastafile)
  find  $dirname -name "$fastaname.harvest" -type f -size -570c -delete 
  find  $dirname -name "$fastaname.harvest.gff3" -type f -size -1c -delete 
  
  sed -i '1,12d' $fastafile.harvest 
    if [ -f "$fastafile.harvest.gff3" ]; then
    cat $fastafile.harvest.gff3 >> $fasta/${process_id}.fna.harvest.combine.gff3
    fi

    name=$(basename $fastafile ".fna")
    sed -i 's/  / /g' $fasta/$name.fna.harvest
    sed -i -r 's/(\s+)?\S+//11' $fasta/$name.fna.harvest
    awk -F "\t" -v chr=$name '{print $0" "'chr'}' $fasta/$name.fna.harvest >>$fasta/$procces_id.all.harvest
 
    sed -i "s/$name/$ID_chromo $NC_chromo/g" $fasta/$procces_id.all.harvest
    $LTRFINDER $fastafile -w 2 -C -D $maxdisltr -d $mindisltr -L $maxlenltr -l $minlenltr -p $matchpairs -M $similarFinder -s $trna >$fastafile.finder
    find  $dirname -name  "$fastaname.finder"  -type f -size -398c -delete ##delet all files .finder less than 398 bytes
    if [ -f "$fastafile.finder" ]; then
        perl $RUN/convert_ltr_finder2.pl  $fastafile.finder $ID_chromo >$fastafile.finder.scn
            if [[ "$ID_chromo" == "0" ]]; then
            head -12 ${fastafile}.finder.scn > $fasta/all.finder.scn.header
            fi
        sed  '1,12d' $fastafile.finder.scn >>$fasta/$procces_id.all.finder.scn
        cat $fastafile.finder >> $fasta/all.fna.finder
    fi
 rm $fastafile