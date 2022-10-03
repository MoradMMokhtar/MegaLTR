#!/usr/bin/env Rscript
library(ggplot2)
library(hrbrthemes)
library(viridis)
require(scales)
options(bitmapType='cairo')

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop()
}

LFILE=args[1]
DFILE=args[2]
PORGRAM=args[3]
plot_functions_path=args[4]

plotscript=paste(plot_functions_path,"plot_functions.r", sep="/")
source(plotscript)

mypath=file.path(LFILE)

### LTR Length
data.len<-read.csv(LFILE,header = F,sep="\t")

data=data.len

valuename="V5"
groupname="V2"
programname=PORGRAM

ppc<-plot_LTR.chart(data,valuename,groupname,"LTR-RT Length (bp)","LTR-RT Types",programname )
ppc
ppb<-plot_LTR.boxplot(data,valuename,groupname,"LTR-RT Length (bp)","LTR-RT Types",programname )
ppb

plot_LTR.save(paste0(mypath,".Length_chart"),ppc)
plot_LTR.save(paste0(mypath,".Length_boxplot"),ppb)


### LTR Length
data.time<-read.csv(DFILE,sep="\t",header = F)

data=data.time
valuename="V3"
groupname="V2"
programname=PORGRAM

ppc<-plot_LTR.chart(data,valuename,groupname,"TimeK (Generation)","LTR-RT Types",programname )
ppc
ppb<-plot_LTR.boxplot(data,valuename,groupname,"TimeK (Generation)","LTR-RT Types",programname )
ppb

plot_LTR.save(paste0(mypath,".TimeK_chart"),ppc)
plot_LTR.save(paste0(mypath,".TimeK_boxplot"),ppb)
