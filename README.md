# # MegaLTR

MegaLTR is a robust online server that identifies intact LTR-RTs in any target genome, and its local standalone version. MegaLTR is freely available at [https://bioinformatics.um6p.ma/MegaLTR](https://bioinformatics.um6p.ma/MegaLTR)


MegaLTR is a pipeline that detects intact LTR-RTs at the whole genome level. The pipeline integrates the structure-based, homology-based and de novo intact LTR-RT identification, classification, annotation and visualization tools such as [LTR_FINEDR](https://github.com/xzhub/LTR_Finder), [LTRharvest](http://genometools.org/pub/binary_distributions/), [LTR_retriever](https://github.com/oushujun/LTR_retriever), [RepeatMasker](http://www.repeatmasker.org/), [CDHIT package](http://weizhongli-lab.org/cd-hit/), [BLAST+ package](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/), [HMMER package](http://hmmer.org/),  [LTRdigest](https://www.zbh.uni-hamburg.de/en/forschung/gi/software/ltrdigest.html), [TEsorter](https://github.com/zhangrengang/TEsorter), l [REANNOTATE](http://www.bioinformatics.org/reannotate/about.html), [ClustalW](https://anaconda.org/bioconda/clustalw), [faidx](https://anaconda.org/bioconda/pyfaidx), [Rscript](https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/Rscript), and [RIdeogram](https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html).


MegaLTR running with three options:

> 1: Intact LTR-RT identification and annotation of internal domains.
> 
> 2: Intact LTR-RT Identification and annotation of internal domains plus determination of insertion time.
> 
> 3: Intact LTR-RT Identification, annotation of internal domains, determination of insertion time, LTR-RT gene-chimera analysis and visualization of gene density and LTR-RTs across chromosomes.

MegaLTR has been tested on Ubuntu 18.04 and 20.04.

 1. Install
 2. Required data
 3. Usage
 4. Run Example
 5. Output files
 6. Output files example
#
## **Install**
The installation require conda. You can install all dependencies for running MegaLTR in a new conda environment using the MegaLTR.yml file. If you do not have conda, please follow [this tutorial](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).

1- Download repository from github 
>`git clone https://github.com/MoradMMokhtar/MegaLTR.git` 

2- Go to the MegaLTR folder 
>`cd MegaLTR ` 

3- Create the MegaLTR environment with all dependencies   
>`conda env create -f MegaLTR.yml`  
  
4- Activate the MegaLTR environment  
>`conda activate MegaLTR`

## Required data

|Data|Option 1|Option 2|Option 3|
|--|--|--|--|
| Genome sequence - Fasta file with chromosomes/scaffolds/contigs sequences | Required |Required  |Required  |
| Genome annotation - GFF/GFF3 file with genome annotations (gene,CDS,mRNA) | Not Required |Not Required  |Required  |
|  |  |  |  |

##  Usage
Go to the MegaLTR folder


    
    bash MegaLTR.sh -A [1 or 2 or 3] -F [Genome in FASTA Format] -G [GFF/GFF3 File]
    
    #Required arguments:
    -A		The analysis type [1 or 2 or 3] 
		    1 (for Intact LTR-RT identification and annotation of internal domains 'This analysis needs FASTA file only') 
		    2 (for Intact LTR-RT Identification and annotation of internal domains plus determination of insertion time 'This analysis needs FASTA file only') 
		    3 (for Intact LTR-RT Identification, annotation of internal domains, determination of insertion time, LTR-RT gene-chimera analysis and visualization of gene density and LTR-RTs across chromosomes 'This analysis needs FASTA and GFF files')     
	-F		Your path to the genome sequence (Fasta file). 
	-G		Your path to the genome annotation (GFF/GFF3 file). Required with argument -A 3 only.

    #Optional arguments:
    -T		tRNA sequence file (Locate the filename from the tRNA folder or provide your own tRNA sequence in FASTA format, default is Arabidopsis_thaliana_trna.fa)."
    -P		Outfileprefix, default is results."
    -l		Min length of 5'&3'LTR, default is 100."
    -L		Max length of 5'&3'LTR, default is 7000."
    -d		Min distance between 5'&3'LTR, default is 1000."
    -D		Max distance between 5'&3'LTR, default is 15000."
    -S		Specify similaritythreshold, default is 85."
    -M		Min length of exact match pair, default is 20."
    -B		TE Database that TEsorter will use it {gydb,rexdb,rexdb-plant,rexdb-metazoa,sine}, default is rexdb."
    -C		Mininum coverage for protein domains in HMMScan output, default is 20."
    -V		Maxinum E-value for protein domains in HMMScan output, default is 0.001."
    -Q		Classifying rule [identity-coverage-length] based on similarity, default is 80-80-80]."
    -E		Hmm-database that TEsorter will use it {gydb,rexdb,rexdb-plant,rexdb-metazoa,sine}, default is rexdb]."
    -R		Neutral mutation rate of the target species (per bp per ya), e.g., rice: 1.3e-8 [0.000000013]; mammal: 2.2e-9 [0.0000000022]; Drosophila: 1.6e-8 [0.000000016], default is 0.000000013."
    -U		The distance upstream LTR retrotransposons, default is 5000."
    -X		The distance downstream LTR retrotransposons, default is 5000."
    -W		Window size to extract gene density from the GFF file, default is 1000000."
    -N		Number of chromosomes specified in FASTA file to visualize density of genes and LTRs, default is 12."
    -v		Print MegaLTR version and exit."
    -t		Indicate how many CPU/threads you want to run MegaLTR, default is 4."
    -h		Print this Help
    
    Default parameters:bash MegaLTR.sh -A 3 -F /path/to/genome_fasta_file  -G /path/to/gff_file -T Arabidopsis_thaliana_trna.fa -P Results -l 100 -L 7000 -d 1000 -D 15000 -S 85 -M 20 -B rexdb -C 20 -V 0.001 -Q 80-80-80 -E rexdb -R 0.000000015 -U 500 -X 5000 -W 1000000 -N 9 -t 6


## Run Example

    bash MegaLTR_Run_Example.sh -A 3 -F NC_003070.9_Arabidopsis_thaliana.fna.gz -G Arabidopsis_thaliana.gff.gz

## Output files

We have collected the main output files in the Collected _Files folder in the main output directory. The results are presented in the form of tables and images as follows:

| # | File name | Description |
|--|--|--|
| 1 | *.fna.pass.list | All LTR-RTs that passed the filtering step |
| 2 | *.fna.nmtf.pass.list | Non-TGCA LTR-RTs that passed the filtering step |
| 3 | *.fna.pass.list.gff3 | GFF3 format for intact LTR-RTs |
| 4 | *.fna.LTRlib.redundant.fa | All LTR-RTs with redundancy in FASTA format |
| 5 | *.fna.LTRlib.fa | All non-redundant LTR-RTs in FASTA format |
| 6 | *.fna.LTR.gff3 | GFF3 format |
| 7 | *.statistics.tsv | LTR family summary |
| 8 | all.finder.scn | The LTR_Finder results |
| 9 | all.harvest | The LTR_harvest results |
| 10 | genes_up_and_down_LTR.tsv | Genes up- and down-stream of LTR-RT elements |
| 11 |LTR_Table_Digest_TEsorter _Time_nongene_and_gene.tsv  | combine the results of LTR-Finder, LTRharvest, LTR-retriever, LTRdigest, TEsorter, insertion time, LTR-RT-gene chimeras, and LTR-RT near genes in one file |
| 12 | LTR_Table_TEsorter_Digest.tsv | Mergeing of LTR_retriever, LTRdigest, and TEsorter results in one file |
| 13 | *.Digest_TEsorter_Time.tsv | Mergeing of LTR_retriever, LTRdigest, TEsorter, and insertion time results in one file |
| 14 | *.fna.harvest.combine.gff3 | The LTR_harvest results in GFF3 format |
| 15 | *.fna.defalse | the false LTR-RTs that does not pass the filtering step |
| 16 | *.fna.masked | All LTRs recognized by RepeatMasker |
| 17 | chromosome.png | visualization of gene density and LTR-RTs across chromosomes PNG format |
| 18 | *.length.ids2.Length_boxplot.png | the boxplot of LTR-RT length for both LTR-RTs superfamilies |
| 19 | *.length.ids2.Length_chart.png | statistical distribution of LTR-RT length for both LTR-RTs superfamilies |
| 20 | *.length.ids2.TimeK_boxplot.png | the boxplot of LTR-RT insertion age for both LTR-RTs superfamilies |
| 21 | *.length.ids2.TimeK_chart.png | statistical distribution of LTR-RT insertion age for both LTR-RTs superfamilies |
| 22 | *.length.ids.Length_boxplot.png | the boxplot of LTR-RT length for each LTR-RTs superfamily  |
| 23 | *.length.ids.Length_chart.png | statistical distribution of LTR-RT length for each LTR-RTs superfamily |
| 24 | *.length.ids.TimeK_boxplot.png | the boxplot of LTR-RT insertion age for each LTR-RTs superfamily |
| 25 | *.length.ids.TimeK_chart.png | statistical distribution of LTR-RT insertion age for each LTR-RTs superfamily |
| 26 | *.PBS.Sequence.fa | All PBS in FASTA format |
| 27 | *.PPT.Sequence.fa | All PPT in FASTA format |
| 28 | LTR-RT_Sequence.fa | All intact LTR-RTs in FASTA format |
 #
 ## Output files example (images)
 
 - Genes and LTR-RTs density across chromosomes

![Genes and LTR-RTs density across chromosomes](https://bioinformatics.um6p.ma/MegaLTR/images/chromosomes1.png)
 - Boxplot of LTR-RTs length for all LTR-RTs

![Boxplot of LTR-RTs length for all LTR-RTs](https://bioinformatics.um6p.ma/MegaLTR/images/new_test_17.length.ids2.Length_boxplot1.png)
 - Statistical distribution of LTR-RTs length for all LTR-RTs
  ![Statistical distribution of LTR-RTs length for all LTR-RTs](https://bioinformatics.um6p.ma/MegaLTR/images/new_test_17.length.ids2.Length_chart1.png)
 - Boxplot of LTR-RTs insertion age for all LTR-RTs
  ![Boxplot of LTR-RTs insertion age for all LTR-RTs](https://bioinformatics.um6p.ma/MegaLTR/images/new_test_17.length.ids2.TimeK_boxplot1.png)
 - Statistical distribution of LTR-RTs insertion age for all LTR-RTs
  ![Statistical distribution of LTR-RTs insertion age for all LTR-RTs](https://bioinformatics.um6p.ma/MegaLTR/images/new_test_17.length.ids2.TimeK_chart1.png)
 -  Boxplot of LTR-RT length for each LTR-RTs class
  ![Boxplot of LTR-RT length for each LTR-RTs class](https://bioinformatics.um6p.ma/MegaLTR/images/new_test_17.length.ids.Length_boxplot1.png)
 - Statistical distribution of LTR-RT length for each LTR-RTs class
![Statistical distribution of LTR-RT length for each LTR-RTs class](https://bioinformatics.um6p.ma/MegaLTR/images/new_test_17.length.ids.Length_chart1.png)
 - Boxplot of LTR-RTs insertion age for each LTR-RTs class
![Boxplot of LTR-RTs insertion age for each LTR-RTs class](https://bioinformatics.um6p.ma/MegaLTR/images/new_test_17.length.ids.TimeK_boxplot1.png)
 - Statistical distribution of LTR-RTs insertion age for each LTR-RTs class
![Statistical distribution of LTR-RTs insertion age for each LTR-RTs class](https://bioinformatics.um6p.ma/MegaLTR/images/new_test_17.length.ids.TimeK_chart1.png)
#
## For more information and help

To report bugs and give us suggestions, you can open an issue here. You may also contact us by e-mail morad.mokhtar@um6p.ma or achraf.elallali@um6p.ma
