import numpy as np
import sys
items = {}

def split_chromo(file_name):
    f = open(file_name, "r")
    global items
    for line in f:
        word = line.split()
        if word[0] in items:
            name = word[0]
            items[name].append(int(word[2]))
        else:
            name = word[0]
            items[name] = []
            items[name].append(int(word[2]))

def get_bins(interval, max):
    data = []
    for x in range(0, max + interval, interval):
        data.append(x)
    return data

def print_result(chr, result, interval, f_export):
    i = 0
    res_len = len(result[0])
    while (i < res_len):
        if (result[0][i] == 0):
            i += 1
            continue
        f_export.write(chr + "\t" + str(result[1][i])+ "\t" + str(result[1][i] + interval) + "\t" + str(result[0][i])+"\n")
        i += 1
    
# Main Function 
def counter(file_name, interval, new_file_name):
    global items
    split_chromo(file_name)
    f_export = open (new_file_name, "w")
    for chr in items:
        items[chr].sort()
        size = len(items[chr])
        print_result(chr, np.histogram(items[chr], bins = get_bins(interval, items[chr][size - 1])), interval, f_export)
    f_export.close()
gene_anno=sys.argv[1]
step=sys.argv[2]
output=sys.argv[3]
#print(gene_anno,step,output)
# Example (file_name, interval, new_file_name)
counter(gene_anno, int(step), output)
