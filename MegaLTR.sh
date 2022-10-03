#!/bin/bash
if [ $# -eq 0 ]
   then
      echo "Parameters -A and -F is required, use [bash MegaLTR.sh -help] for more detalis"
      exit
fi
tRNAdb=$(pwd)/bin/tRNA
LTR_FINDER=$(pwd)/bin/LTR_FINDER.x86_64-1.0.7/ltr_finder
chmod 775 $LTR_FINDER #Give execulting permissions for ltr_finder
LTR_FINDER_convert=$(pwd)/bin/RUN/convert_ltr_finder2.pl
LTRretriever=$(pwd)/bin/LTR_retriever/LTR_retriever
chmod 775 $LTRretriever #Give execulting permissions for LTRretriever
RUN=$(pwd)/bin/RUN
eval "$(conda shell.bash hook)"
############################################################
# Set variables                                            #
############################################################
Analysistype=3
Fastafilepath=""
gffpath=""
trna=Arabidopsis_thaliana_trna.fa
outfileprefix=results
minlenltr=100
maxlenltr=7000
mindisltr=1000
maxdisltr=15000
similar=85
matchpairs=20
TEsorterdb=rexdb
TEsortercov=20
TEsortereval=0.001
TEsorterrule=80-80-80
TEsorterhmm=rexdb
RateOfEvolution=0.000000015
up=500
down=5000
density1=1000000
numberofchromosom=9
threads=4
############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "MegaLTR is a pipeline designed for high-throughput identification and classification of LTR-retrotransposons, insertion age estimation, detection of LTR-gene chimeras, and determination of nearby genes."
   echo
   echo "     MegaLTR v1.07."
   echo
   echo "     Options:"
   echo "     -h     Print this Help."
   echo "     -A     The analysis type [1 or 2 or 3] 1 (for Intact LTR-RT identification and annotation of internal domains 'This analysis needs FASTA file only') 2 (for Intact LTR-RT Identification and annotation of internal domains plus determination of insertion time 'This analysis needs FASTA file only') 3 (for Intact LTR-RT Identification, annotation of internal domains, determination of insertion time, LTR-RT gene-chimera analysis and visualization of gene density and LTR-RTs across chromosomes 'This analysis needs FASTA and GFF files') , default is 3"
   echo "     -F     Your path to Fasta file."
   echo "     -G     Your path to GFF file."
   echo "     -T     tRNA sequence file (Locate the filename from the tRNA folder or provide your own tRNA sequence in FASTA format, default is Arabidopsis_thaliana_trna.fa)."
   echo "     -P     Outfileprefix, default is results."
   echo "     -l     Min length of 5'&3'LTR, default is 100."
   echo "     -L     Max length of 5'&3'LTR, default is 7000."
   echo "     -d     Min distance between 5'&3'LTR, default is 1000."
   echo "     -D     Max distance between 5'&3'LTR, default is 15000."
   echo "     -S     Specify similaritythreshold, default is 85."
   echo "     -M     Min length of exact match pair, default is 20."
   echo "     -B     TE Database that TEsorter will use it {gydb,rexdb,rexdb-plant,rexdb-metazoa,sine}, default is rexdb."
   echo "     -C     mininum coverage for protein domains in HMMScan output, default is 20."
   echo "     -V     maxinum E-value for protein domains in HMMScan output, default is 0.001."
   echo "     -Q     classifying rule [identity-coverage-length] based on similarity, default is 80-80-80]."
   echo "     -E     hmm-database that TEsorter will use it {gydb,rexdb,rexdb-plant,rexdb-metazoa,sine}, default is rexdb]."
   echo "     -R     neutral mutation rate of the target species (per bp per ya), e.g., rice: 1.3e-8 [0.000000013] (Ma and Bennetzen 2004); mammal: 2.2e-9 [0.0000000022] (S. Kumar 2002); Drosophila: 1.6e-8 [0.000000016](Bowen and McDonald 2001), default is 0.000000013."
   echo "     -U     the distance upstream LTR retrotransposons, default is 5000."
   echo "     -X     the distance downstream LTR retrotransposons, default is 5000."
   echo "     -W     window size to extract gene density from the GFF file, default is 1000000."
   echo "     -N     Number of chromosomes specified in FASTA file to visualize density of genes and LTRs, default is 12."
   echo "     -v     Print MegaLTR version and exit."
   echo "     -t     Indicate how many CPU/threads you want to run MegaLTR, default is 4."
   echo
   echo "Default parameters:bash MegaLTR.sh -A $Analysistype -F Your_path_to_FASTA_file -G Your_path_to_GFF_file -T $trna -P $outfileprefix -l $minlenltr -L $maxlenltr -d $mindisltr -D $maxdisltr -S $similar -M $matchpairs -B $TEsorterdb -C $TEsortercov -V $TEsortereval -Q $TEsorterrule -E $TEsorterhmm -R $RateOfEvolution -U $up -X $down -W $density1 -N $numberofchromosom -t $threads"
   exit
}
version()
{
   echo "MegaLTR v1.06."
   echo
   exit
}
############################################################
# Process the input options.                               #
############################################################
check=0
while getopts h:A:F:G:T:P:l:L:d:D:S:M:B:C:V:Q:E:R:U:X:W:N:v:t: options
do
   case $options in
      h) Help;;
      v) version;;
      A) Analysistype=$OPTARG;((check=check+1));;
      F) Fastafilepath=$OPTARG;((check=check+1));;
      G) gffpath=$OPTARG;;
      T) trna=$OPTARG;;
      P) outfileprefix=$OPTARG;;
      l) minlenltr=$OPTARG;;
      L) maxlenltr=$OPTARG;;
      d) mindisltr=$OPTARG;;
      D) maxdisltr=$OPTARG;;
      S) similar=$OPTARG;;
      M) matchpairs=$OPTARG;;
      B) TEsorterdb=$OPTARG;;
      C) TEsortercov=$OPTARG;;
      V) TEsortereval=$OPTARG;;
      Q) TEsorterrule=$OPTARG;;
      E) TEsorterhmm=$OPTARG;;
      R) RateOfEvolution=$OPTARG;;
      U) up=$OPTARG;;
      X) down=$OPTARG;;
      W) density1=$OPTARG;;
      N) numberofchromosom=$OPTARG;;
      t) threads=$OPTARG;;
     *) echo "Error: Invalid option" ;;
   esac
done

########################################################
division=100
similarFinder=$(perl -e "print $similar/$division") ###convert the similarity to LTR_FINDER format
process_id=$outfileprefix
userpath= mkdir -p $(pwd)/$process_id
userpath=$(pwd)/$process_id
LAI= mkdir -p $userpath/LAI
LAI=$userpath/LAI
ltrdigest= mkdir -p $userpath/LTRdigest
ltrdigest=$userpath/LTRdigest
mkdir -p $userpath/TEsorter
TEsorter=$userpath/TEsorter
mkdir -p $userpath/Others
Others=$userpath/Others
mkdir -p $userpath/Time
time=$userpath/Time
mkdir -p $userpath/R_plots
Rplots=$userpath/R_plots
mkdir -p $userpath/inside_genes
inside_genes=$userpath/inside_genes
mkdir -p $userpath/near_genes
near_genes=$userpath/near_genes
mkdir -p $userpath/density
densitypath=$userpath/density
process_id=$outfileprefix
mkdir -p $userpath/Collected_Files
Collected_Files=$userpath/Collected_Files
mkdir -p $userpath/FASTA
FASTA=$userpath/FASTA
mkdir -p $userpath/LTR_FINDER
FINDER=$userpath/LTR_FINDER
mkdir -p $userpath/LTR_HARVEST
LTRHARVEST=$userpath/LTR_HARVEST

###############################################################
#### MegaLTR process the start.                               #
###############################################################
now="$(date)" 
echo
printf "\t#############################################\n\t##############  MegaLTR v1.06  ##############\n\t#############################################\n\n\tContributors: Morad M Mokhtar, Achraf El Allali\n\n"
printf "\t$now \t Start time %s\n"  ### print current date
echo
printf "\tParameters: -A $Analysistype -F $Fastafilepath -G $gffpath -T $trna -P $outfileprefix -l $minlenltr -L $maxlenltr -d $mindisltr -D $maxdisltr -S $similar -M $matchpairs -B $TEsorterdb -C $TEsortercov -V $TEsortereval -Q $TEsorterrule -E $TEsorterhmm -R $RateOfEvolution -U $up -X $down -W $density1 -N $numberofchromosom -t $threads\n\n"

if [ $check -ne 2 ]
   then
      printf "\tParameters -A and -F is required, use [bash MegaLTR.sh -help] for more detalis\n"
      exit
fi
      conda activate MegaLTR
if [[ $Fastafilepath =~ \.gz$ ]]; then
   cd $userpath
      echo
      printf "\tCheck the FASTA File format.\n"
      gunzip -c "$Fastafilepath" >$userpath/$process_id.fna
      perl $RUN/checkfasta.pl $userpath/$process_id.fna ### Check the FASTA File format
   elif [[ $Fastafilepath =~ \.zip$ ]]; then
      echo
      printf "\tCheck the FASTA File format.\n"
      gunzip -c "$Fastafilepath" >$userpath/$process_id.fna
      perl $RUN/checkfasta.pl $userpath/$process_id.fna ### Check the FASTA File format
   elif ([ $(stat -c%s "$Fastafilepath") -gt 500 ]); then
      printf "\tCheck the FASTA File format.\n"
      perl $RUN/checkfasta.pl $Fastafilepath ### Check the FASTA File format
      cp $Fastafilepath $userpath/$process_id.fna #### copy fasta file tRNA files
   else
      printf "\n\tCheck the FASTA File format.\n"
      printf "\n\tPlease use -F to provide a valid FASTA file [-F FASTA_file_path].\n"
      exit
fi

if ([ $Analysistype -eq 3 ] && [[ $gffpath == "" ]]); then
      printf "\n\tPlease use -G to provide a valid GFF file [-G GFF_file_path].\n"
   exit
fi
if ([[ $Analysistype -eq 3 ]] && [ $(stat -c%s "$gffpath") -lt 500 ]); then
      printf "\n\tPlease use -G to provide a valid GFF file [-G GFF_file_path].\n"
   exit
fi

if ([ $Analysistype -eq 3 ] && [[ $gffpath =~ \.gz$ ]]); ### mode 3
   then
      gunzip -c "$gffpath" >$userpath/$process_id.gff
   elif ([ $Analysistype -eq 3 ] && [[ $gffpath =~ \.zip$ ]]); then
      gunzip -c "$gffpath" >$userpath/$process_id.gff
   elif ([[ $Analysistype -eq 3 ]] && [ $(stat -c%s "$gffpath") -gt 500 ]); then
      cp $gffpath $userpath/$process_id.gff #### copy gff file
fi
if ([ $(stat -c%s "$userpath/$process_id.fna") -gt 500 ]);
   then
      ## {  ############## copy and split fasta #############
         conda activate MegaLTR
            cp $tRNAdb/$trna $userpath/$trna  #### copy tRNA file 
            sed -i 's/ .*//' $userpath/$process_id.fna  ### remove unnecessary details from Fasta sequence headers
            sed -i 's/\t.*//g' $userpath/$trna  ### remove unnecessary details from tRNA sequence headers
            cd $FASTA
            now2="$(date)"
            echo
            printf "\n\t$now2 \tLTR_FINDER Started %s\n"
            faidx --split-files $userpath/$process_id.fna
      ## }
      ## #{    ################## LTR_Finder ########### 
            python3 $RUN/LTR_FINDER_threads.py $FASTA $threads $LTR_FINDER $maxdisltr $mindisltr $maxlenltr $minlenltr $matchpairs $similarFinder $userpath $trna $LAI $FASTA
            find . -name "*.finder" -type f -size -398c -delete ##delet all files .finder less than 398 bytes
            cat $FASTA/*.fna.finder >$FASTA/all.fna.finder
            perl $RUN/convert_ltr_finder2.pl $FASTA/all.fna.finder > $FASTA/all.finder.scn
            sed -i -r 's/(\s+)?\S+//11' $FASTA/all.finder.scn
            sed '1,12d' $FASTA/all.finder.scn > $FASTA/all.finder.scn2
            cp $FASTA/all.finder.scn $Collected_Files
            sed -i '1,7d' $FASTA/*.fna.finder
            sed -i '$ d' $FASTA/*.fna.finder
            sed -i 's/No LTR Retrotransposons Found//g' $FASTA/*.fna.finder
            cat $FASTA/*.fna.finder >$FINDER/$process_id.fna.LTR_finder.1
            echo "index	SeqID	Location	LTR len	Inserted element len	TSR	PBS	PPT	RT	IN (core)	IN (c-term)	RH	Strand	Score	Sharpness	Similarity" >$FINDER/$process_id.fna.LTR_finder.tsv
            awk 'NF > 0' $FINDER/$process_id.fna.LTR_finder.1 >>$FINDER/$process_id.fna.LTR_finder.tsv
            cut -f2- $FINDER/$process_id.fna.LTR_finder.tsv >$FINDER/$process_id.fna.LTR_finder.tsv2
            sed -E 's/_sub\S+\t/\t/g' $FINDER/$process_id.fna.LTR_finder.tsv2 >$FINDER/$process_id.fna.LTR_finder.tsv
            now3="$(date)"
            echo
            printf "\t$now3 \tLTR_FINDER done %s\n"
      ## #}
#   # # {  ################ Run LTR_HARVEST #############
            now4="$(date)"
            echo
            printf "\t$now4 \tLTR_HARVEST Started %s\n"
		      for fst in  $FASTA/*.fna 
		      do
   			   gt suffixerator -db $fst -indexname $fst -tis -suf -lcp -des -ssp -sds -dna
		      done
		      python3 $RUN/LTR_HARVEST_threads.py $FASTA $threads $maxdisltr $mindisltr $maxlenltr $minlenltr $matchpairs $similar $LTRHARVEST
            rm $FASTA/*.des $FASTA/*.esq $FASTA/*.lcp $FASTA/*.llv $FASTA/*.md5 $FASTA/*.prj $FASTA/*.sds $FASTA/*.suf
            find . -name "*.harvest" -type f -size -570c -delete ##delet all files .harvest less than 570 bytes
            find . -name "*.harvest.gff3" -type f -size -1c -delete ##delet all files .harvest.gff3 less than 1 bytes
            sed -i '1,12d' $FASTA/*.harvest
            for i in $FASTA/*.harvest
            do
               name=$(basename $i ".fna.harvest")
               sed -i 's/  / /g' $FASTA/$name.fna.harvest
               sed -i -r 's/(\s+)?\S+//11' $FASTA/$name.fna.harvest
               awk -F "\t" -v chr=$name '{print $0" "'chr'}' $i >$i.2
            done
            cat $FASTA/*.harvest.2 >$FASTA/all.harvest.3
            printf "# args=-index $FASTA/$process_id.fna -maxdistltr $maxdisltr -mindistltr $mindisltr -maxlenltr $maxlenltr -minlenltr $minlenltr -similar $similar -seqids yes -gff3 $FASTA/$process_id.fna.harvest.gff3\n# predictions are reported in the following way\n# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr \n# where:\n# s = starting position\n# e = ending position\n# l = length\n# ret = LTR-retrotransposon\n# lLTR = left LTR\n# rLTR = right LTR\n# sim = similarity\n# seq-nr = sequence number" >$FASTA/all.harvest
            cat $FASTA/*.harvest.3 >>$FASTA/all.harvest
            cat $FASTA/all.harvest $FASTA/all.finder.scn >$FASTA/all.harvest.finder.combine
            awk -F " " '{print $11}' $FASTA/all.harvest.3 > $FASTA/ids.list
            awk -F " " '{print $11}' $FASTA/all.finder.scn2 >> $FASTA/ids.list
            awk '!a[$0]++' $FASTA/ids.list >$FASTA/ids.uniq  
            awk '{print $0 "\t" NR-1}' $FASTA/ids.uniq >$FASTA/ids.uniq_with_numbers
            while read a b
            do
               sed -i "s/$a/$b $a/g" $FASTA/all.harvest.finder.combine
            done < $FASTA/ids.uniq_with_numbers
            cp $FASTA/all.harvest.finder.combine $LAI
            cp $FASTA/all.harvest $Collected_Files
            cat $FASTA/*.harvest.gff3 >$FASTA/$process_id.fna.harvest.combine.gff3
            cp $FASTA/$process_id.fna.harvest.combine.gff3 $Collected_Files
            now5="$(date)"
            echo
            printf "\t$now5 \tLTR_HARVEST Done %s\n"
#   # # }
#  # # {     #######LTR_retriever ##########

         cd $LAI
         now6="$(date)"
         echo
         printf "\t$now6 \tLTR_retriever Started %s\n"
         bash $RUN/newLTR_retriever.sh  $userpath/$process_id.fna $userpath $LAI $minlenltr $process_id $Collected_Files $threads $LTRretriever
         now7="$(date)"
         echo
         printf "\t$now7 \tLTR_retriever Done %s\n"

#   # # }
  # # {     #######ltrdigest ###########

conda activate MegaLTR
      now8="$(date)"
      echo
      printf "\t$now8 \tLTRdigest Started %s\n"
      bash $RUN/ltrdigest.sh $process_id $userpath $trna $LAI $ltrdigest
      now9="$(date)"
      echo
      printf "\t$now9 \tLTRdigest Done %s\n"
  # # }
  # # {     #######TEsorter ###########
      now10="$(date)"
      echo
      printf "\t$now10 \tTEsorter Started %s\n"
      cd $TEsorter
      bash $RUN/TEsorter.sh $TEsorter $ltrdigest/"$process_id"_complete.fas $TEsorterdb $TEsortercov $TEsortereval $TEsorterrule $TEsorterhmm $RUN $process_id $Collected_Files $userpath $threads
      now11="$(date)"
      echo
      printf "\t$now11 \tTEsorter Done %s\n"
      now12="$(date)"
      echo
      printf "\t$now12 \tFiltering TEsorter results Started %s\n"
      awk -F '\t' '$2 == "LTR"' $TEsorter/"$process_id"_complete.fas.$TEsorterhmm.cls.tsv >$Others/$process_id.LTR.tsv  ##### search inside a specific column
      awk -F '\t' '$2 != "LTR"' $TEsorter/"$process_id"_complete.fas.$TEsorterhmm.cls.tsv >$Others/$process_id.others.tsv  ##### search inside a specific column
      perl $RUN/print.pl $ltrdigest/"$process_id"_tabout.csv >$Others/$process_id.tabout.tsv ###### $Others/$"$process_id".tabout.tsv header ########## element	id element start	element end	element length	sequence	lLTR start	lLTR end	lLTR length	rLTR start	rLTR end	rLTR length	lTSD start	lTSD end	lTSD motif	rTSD start	rTSD end	rTSD motif	PPT start	PPT end	PPT motif	PPT strand	PPT offset	PBS start	PBS end	PBS strand	tRNA	tRNA motif	PBS offset	tRNA offset	PBS/tRNA edist
      now13="$(date)"
      echo
      printf "\t$now13 \tMergeing of LTR_retriever, LTRdigest, and TEsorter results %s\n"
      perl $RUN/TEsorter_Digest.pl  $Others/$process_id.tabout.tsv $Others/$process_id.LTR.tsv >$Others/LTR_Table_TEsorter_Digest.tsv
      cp $Others/LTR_Table_TEsorter_Digest.tsv $Collected_Files ### LTR_retriever, LTRdigest, and TEsorter results in one file 
  # # }   
  # # {

      if ([ $Analysistype -eq 2 ] || [ $Analysistype -eq 3 ]) ### mode 2
         then
conda activate MegaLTR
            ########LTR insertion time ###########
            now14="$(date)"
            echo
            printf "\t$now14 \tCalculation of the LTR-RT insertion time %s\n"
            python3 $RUN/get_region.py $userpath/$process_id.fna $Others/LTR_Table_TEsorter_Digest.tsv $time
            echo "" > $Others/$process_id.time.txt
            echo "LTRNAME	K_Kimura	K_Ksd	K_TajimaNei	K_TNsd	timeK	timeKsd	numComparedSites	transitions	numComparedSites	transversions	numComparedSites	timeTN	transitions_numComparedSites	transversions_numComparedSites" >>$Others/$process_id.time.txt
               for fst in  $time/*.fasta 
               do
                  clustalw -infile="$fst" >>$Others/$process_id.clustalw.txt  ### Alignment of the LTR using clustalw
                  fname=$(basename $fst ".fasta")
                  ltrk=$(perl $RUN/estimate_K.pl $time/$fname".aln" $RateOfEvolution) ## estimate LTR-RT insertion time using Kimura and Tajima&Nei methods (this script is part of REannotate program)
                  echo -e $fname "\t""$ltrk"  >> $Others/$process_id.time.txt
               done
            perl $RUN/TEsorterandtable_time.pl  $Others/LTR_Table_TEsorter_Digest.tsv $Others/$process_id.time.txt >$Others/$process_id.Digest_TEsorter_Time.tsv  ############# $Others/$process_id.Digest_TEsorter_Time (header) org_name	acc	element_start	element_end	element_length			strand	lTSD_start	lTSD_end	lLTR_start	lLTR_end	rLTR_start	rLTR_end	rTSD_start	rTSD_end	EDTA_infoo	PPT start	PPT end	PPT motif	PPT offset	PBS start	PBS end	tRNA	tRNA motif	PBS offset	tRNA offset	PBS/tRNA	edist Order	Superfamily	Clade	Complete	Strand	Domains	LTRNAME	K_Kimura	K_Ksd	K_TajimaNei	K_TNsd	timeK	timeKsd	numComparedSites	transitions	numComparedSites	transversions	numComparedSites	timeTN	transitions_numComparedSites	transversions_numComparedSites
            cp $Others/$process_id.Digest_TEsorter_Time.tsv $Collected_Files ### LTR_retriever, LTRdigest, TEsorter, and insertion time results in one file
            ###### preparing the files for R plots #
            now15="$(date)"
            echo
            printf "\t$now15 \tPreparing R plots files %s\n"
            awk -F "\t" '{ print  $1"\t"$33"\t"$42 }' $Others/$process_id.Digest_TEsorter_Time.tsv > $Others/$process_id.time.ids1
            awk -F '\t' '$2 == "Copia"' $Others/$process_id.time.ids1 >>$Rplots/$process_id.time.ids 
            awk -F '\t' '$2 == "Gypsy"' $Others/$process_id.time.ids1 >>$Rplots/$process_id.time.ids 
            awk -F "\t" '{ print  $2"\t"$33"\t"$3"\t"$4"\t"$5 }' $Others/$process_id.Digest_TEsorter_Time.tsv > $Others/$process_id.length.ids1
            awk -F '\t' '$2 == "Copia"' $Others/$process_id.length.ids1 >>$Rplots/$process_id.length.ids
            awk -F '\t' '$2 == "Gypsy"' $Others/$process_id.length.ids1 >>$Rplots/$process_id.length.ids
            awk -F "\t" '{ print  $1"\t""Copia and Gypsy""\t"$3 }' $Rplots/$process_id.time.ids > $Rplots/$process_id.time.ids2
            awk -F "\t" '{ print  $1"\t""Copia and Gypsy""\t"$2"\t"$3"\t"$4"\t"$5 }' $Rplots/$process_id.length.ids >$Rplots/$process_id.length.ids2
            ##module load R
            Rscript $RUN/chart_plot.r $Rplots/$process_id.length.ids $Rplots/$process_id.time.ids  MegaLTR..$process_id $RUN 2>/dev/null ### Both LTR-RT insertion time and LTR-RT length destributions plots based on the LTR-RT superfamilies
            Rscript $RUN/chart_plot.r $Rplots/$process_id.length.ids2 $Rplots/$process_id.time.ids2  MegaLTR..$process_id  $RUN 2>/dev/null ### Both LTR-RT insertion time and LTR-RT length destributions plots for the whole genome
            cp $Rplots/*.png $Collected_Files ## Both LTR-RT insertion time and LTR-RT length destributions plots (.png)
            now16="$(date)"
            echo
            printf "\t$now16 \tLTR-RT insertion time plots done %s\n"
         ##exit
      fi
  # # }
  # # {   

      if [ $Analysistype -eq 3 ] ### mode 3      # #####  TE-gene chimeras #########
         then
conda activate MegaLTR
            now17="$(date)"
            echo
            printf "\t$now17 \tLTR-RT-gene chimeras started %s\n"
            grep -P "\tgene\t" $gffpath > $userpath/$process_id.grep.gene ## retrieve gene start and end from GFF file
            grep -P "\tpseudogene\t" $gffpath > $userpath/$process_id.grep.pseudogene ## retrieve pseudogene start and end from GFF file
            cat $userpath/$process_id.grep.gene  $userpath/$process_id.grep.pseudogene > $userpath/$process_id.gene_pseudogene.gff  ## combine gene and pseudogene in one file
            perl $RUN/get-TE-within-gene.pl $Others/$process_id.Digest_TEsorter_Time.tsv  $userpath/$process_id.gene_pseudogene.gff $process_id $inside_genes > $inside_genes/$process_id.LTR_inside_genes.table  ## determine the LTR-RT located inside the gene start and end
            awk -F "\t" '{ print $1 }' $inside_genes/$process_id.LTR_inside_genes.table >$inside_genes/$process_id.LTR_inside_genes.table.ids
            awk -F "\t" '{ print $1 }' $Others/$process_id.Digest_TEsorter_Time.tsv >$inside_genes/$process_id.LTR_Table_Digest_TEsorter_Time.ids
            grep -F -x -v -f $inside_genes/$process_id.LTR_inside_genes.table.ids $inside_genes/$process_id.LTR_Table_Digest_TEsorter_Time.ids >$inside_genes/LTR_Table_Digest_TEsorter_Time_process   ###### print found in file 2 an not found in file 1 
            grep -f $inside_genes/LTR_Table_Digest_TEsorter_Time_process $Others/$process_id.Digest_TEsorter_Time.tsv > $inside_genes/LTR_Table_Digest_TEsorter_Time_nongene1
            perl $RUN/No.pl $inside_genes/LTR_Table_Digest_TEsorter_Time_nongene1 >$inside_genes/LTR_Table_Digest_TEsorter_Time_nongene  ## determine which LTR-RT located outside the gene start and end
            cat $inside_genes/$process_id.LTR_inside_genes.table $inside_genes/LTR_Table_Digest_TEsorter_Time_nongene >$inside_genes/LTR_Table_Digest_TEsorter_Time_nongene_and_gene.tsv ## combine the LTR-RT located insid and outside the gene start and end in one file
            cp $inside_genes/LTR_Table_Digest_TEsorter_Time_nongene_and_gene.tsv $Collected_Files
            now18="$(date)"
            echo
            printf "\t$now18 \tLTR-RT-gene chimeras done %s\n"
            #####*** TE near genes***
            now19="$(date)"
            echo
            printf "\t$now19 \tLTR-RT near genes started %s\n"         
            sed  's/\t?\t/\t+\t/g' $userpath/$process_id.gene_pseudogene.gff >$near_genes/$process_id.gene_pseudogene.gff2
            sed  's/\t?\t/\t+\t/g' $Others/$process_id.Digest_TEsorter_Time.tsv >$near_genes/$process_id.Digest_TEsorter_Time2
            awk -F '\t' '$36=="+"'  $near_genes/$process_id.Digest_TEsorter_Time2 >$near_genes/$process_id.Digest_TEsorter_Time5+
            awk -F '\t' '$36=="-"'  $near_genes/$process_id.Digest_TEsorter_Time2 >$near_genes/$process_id.Digest_TEsorter_Time5-
            perl $RUN/get-TE-near-gene-minuse1k.pl $near_genes/$process_id.Digest_TEsorter_Time5+ $near_genes/$process_id.gene_pseudogene.gff2 $up $down >$near_genes/up_genes+
            perl $RUN/get-TE-near-gene-pluse-1k.pl $near_genes/$process_id.Digest_TEsorter_Time5+ $near_genes/$process_id.gene_pseudogene.gff2 $up $down >$near_genes/down_genes+
            perl $RUN/get-TE-near-gene-minuse-1k.pl $near_genes/$process_id.Digest_TEsorter_Time5- $near_genes/$process_id.gene_pseudogene.gff2 $up $down >$near_genes/down_genes-
            perl $RUN/get-TE-near-gene-pluse+1k.pl $near_genes/$process_id.Digest_TEsorter_Time5- $near_genes/$process_id.gene_pseudogene.gff2 $up $down >$near_genes/up_genes-
            cat $near_genes/up_genes+ $near_genes/down_genes+ $near_genes/down_genes- $near_genes/up_genes- >$near_genes/genes_up_and_down_LTR.tsv
            cp $near_genes/genes_up_and_down_LTR.tsv $Collected_Files
            now20="$(date)"
            echo
            printf "\t$now20 \tLTR-RT near genes done %s\n"
            ################### For chromosome LTR-RT distribution ###########
            now20="$(date)"
            echo
            printf "\t$now20 \tVisualization of gene density and LTR-RTs across chromosomes %s\n"
            awk -F "\t" '{ print  $33"\t"$33"\t"$2"\t"$3"\t"$4"\t"$33 }' $Others/$process_id.Digest_TEsorter_Time.tsv > $densitypath/$process_id.figure.distrbution1
            awk -F '\t' 'BEGIN { OFS=FS } { gsub("Copia", "box", $2); print }' $densitypath/$process_id.figure.distrbution1 > $densitypath/$process_id.figure.distrbution2
            awk -F '\t' 'BEGIN { OFS=FS } { gsub("Copia", "6a3d9a", $6); print }' $densitypath/$process_id.figure.distrbution2 > $densitypath/$process_id.figure.distrbution3
            awk -F '\t' 'BEGIN { OFS=FS } { gsub("Gypsy", "triangle", $2); print }' $densitypath/$process_id.figure.distrbution3 > $densitypath/$process_id.figure.distrbution4
            awk -F '\t' 'BEGIN { OFS=FS } { gsub("Gypsy", "ff7f00", $6); print }' $densitypath/$process_id.figure.distrbution4 > $densitypath/$process_id.figure.distrbution5
            grep "##sequence-region" $gffpath > $densitypath/"$process_id"_karyotype
            sed -i 's/##sequence-region //g' $densitypath/"$process_id"_karyotype
            sed -i 's/ /\t/g' $densitypath/"$process_id"_karyotype
            awk -F "\t" '{print $1"\t"$4"\t"$5"\t""1"}' $userpath/$process_id.gene_pseudogene.gff >$densitypath/gene_anno
            python3 $RUN/counter.py $densitypath/gene_anno $density1 $densitypath/gene_anno_counter
            ulimit -s 15947772
            cd $densitypath
               if [ $numberofchromosom -le 8 ] ### the numbers of chromosomes less than or equal 8 to justfay the image
                  then  
                  sort  -k3,3nr $densitypath/"$process_id"_karyotype >$densitypath/"$process_id"_karyotype_sort
                  echo "Chr	Start	End" >$densitypath/"$process_id"_karyotype_sort_head
                  head -8 $densitypath/"$process_id"_karyotype_sort >>$densitypath/"$process_id"_karyotype_sort_head
                  awk -F "\t" '{print $1}' $densitypath/"$process_id"_karyotype_sort_head >$densitypath/"$process_id"_karyotype_sort_head.ids
                  echo "Chr	Start	End	Value" >$densitypath/gene_anno_counter_ok
                  grep -f $densitypath/"$process_id"_karyotype_sort_head.ids $densitypath/gene_anno_counter >>$densitypath/gene_anno_counter_ok
                  echo "Type	Shape	Chr	Start	End	color" >$densitypath/figure_distrbution5_ok
                  grep -f $densitypath/"$process_id"_karyotype_sort_head.ids $densitypath/$process_id.figure.distrbution5 >>$densitypath/figure_distrbution5_ok
                  Rscript $RUN/density_width.r $densitypath/"$process_id"_karyotype_sort_head $densitypath/gene_anno_counter_ok $densitypath/figure_distrbution5_ok
                  cp $densitypath/chromosome.svg $Collected_Files
                  cp $densitypath/chromosome.png $Collected_Files
               fi
               if [ $numberofchromosom -ge 8 ] ### the numbers of chromosomes more than or equal 8 to justfay the image
                  then
                  sort  -k3,3nr $densitypath/"$process_id"_karyotype >$densitypath/"$process_id"_karyotype_sort
                  echo "Chr	Start	End" >$densitypath/"$process_id"_karyotype_sort_head
                  head -20 $densitypath/"$process_id"_karyotype_sort >>$densitypath/"$process_id"_karyotype_sort_head
                  awk -F "\t" '{print $1}' $densitypath/"$process_id"_karyotype_sort_head >$densitypath/"$process_id"_karyotype_sort_head.ids
                  echo "Chr	Start	End	Value" >$densitypath/gene_anno_counter_ok
                  grep -f $densitypath/"$process_id"_karyotype_sort_head.ids $densitypath/gene_anno_counter >>$densitypath/gene_anno_counter_ok
                  echo "Type	Shape	Chr	Start	End	color" >$densitypath/figure_distrbution5_ok
                  grep -f $densitypath/"$process_id"_karyotype_sort_head.ids $densitypath/$process_id.figure.distrbution5 >>$densitypath/figure_distrbution5_ok
                  Rscript $RUN/density.r $densitypath/"$process_id"_karyotype_sort_head $densitypath/gene_anno_counter_ok $densitypath/figure_distrbution5_ok
                  cp $densitypath/chromosome.svg $Collected_Files
                  cp $densitypath/chromosome.png $Collected_Files            
               fi

         ##exit
      fi
  # #}
      rm $userpath/$process_id.fna
      rm $FASTA/*.fna
      now21="$(date)"
      echo
      printf "\t$now21 \tMegaLTR Done, The results saved in ($Collected_Files) %s\n"
   exit 1
fi