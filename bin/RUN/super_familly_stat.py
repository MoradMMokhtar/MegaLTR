#!/usr/bin/env python
# coding: utf-8
import sys,os
infile=sys.argv[1]
outfile=sys.argv[2]
data={}
with open (infile,"r") as file:
    for line in file.readlines():
        splitedline=line.split("\t")
        if splitedline[1] not in data.keys():
            data[splitedline[1]]=[]
        data[splitedline[1]].append(int(splitedline[4].rstrip()))
with open(outfile,'w') as outfile:
    outfile.write(f'LTR-RT class\tNumber of elements\tTotal size bp\n')
    for key in data:
        outfile.write(f'{key}\t{len(data[key])}\t{sum(data[key])}\n')
