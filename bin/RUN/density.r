library("RIdeogram")
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop()
}

karyotype=args[1]
gene_density=args[2]
LTR=args[3]

chr_karyotype <- read.table(karyotype, sep = "\t", header = T, stringsAsFactors = F)
gene_density <- read.table(gene_density, sep = "\t", header = T, stringsAsFactors = F)
Random_RNAs_500 <- read.table(LTR, sep = "\t", header = T, stringsAsFactors = F)

ideogram(karyotype = chr_karyotype, overlaid = gene_density, label = Random_RNAs_500, label_type = "marker")
convertSVG("chromosome.svg", device = "png")
