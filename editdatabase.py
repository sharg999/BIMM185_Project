# TCGA file parser


__author__ = 'Sharong'
#import numpy as np
#from scipy.stats import ttest_1samp, wilcoxon, ttest_ind, mannwhitneyu
import csv
import glob

import os
import fnmatch

#PATH = "C:\Users\Sharong\Desktop\Backed up Data 12.14.2012\Desktop\UCSD\2013-2014 academic year\BIMM 185\DB\BRCA\miRNASeq\BCGSC__IlluminaHiSeq_miRNASeq\Level_3"


def getFiles():
    print("Location of files: ")
    directory = input()

    path = r"%s" % directory
    #dictionary for normal and tumor samples
    allcountsN = dict()
    allcountsT = dict()

    #counts how many miRNA data I have per each
    howmanyN = dict()
    howmanyT = dict()

    #01 at a specific location means tumor sample
    tpat = "01"

    for file in os.listdir(path):
        #parse files for tumor samples
        if fnmatch.fnmatch(file, '*.mirna.quantification.*') and (tpat==file[13:15]):
            #print(file[13:15])
            f = open('C:\\Users\\Sharong\\Desktop\\Level_3\\%s' % file, 'r')
            lines = f.readlines()
            #print(lines)
            for l in lines:
                line = l.split('\t')
                line= list(line)
                #if it's the file title just ignore
                if line[0].startswith("miRNA_ID"):
                    continue
                #if it's the first line of data, ignore the initial string
                elif line[0].startswith("crossed-mapped"):
                    temp = line[0][14:]
                    allcountsT[temp]= float(line[2])
                    howmanyT[temp] = 1
                    #first position in allvals list will be the miRNA_ID
                    #allvalsT.append(temp)
                    #allvalsT.append(float(line[2]))
                #otherwise, just add to dictionary and keep track of sample counts
                elif line[0] in allcountsT:
                    allcountsT[line[0][1:]]+= float(line[2])
                    howmanyT[line[0][1:]]+=1
                    #allvalsT.append(float(line[2]))
                else:
                    allcountsT[line[0][1:]]= float(line[2])
                    howmanyT[line[0][1:]] = 1
                    #allvalsT.append(float(line[2]))

        #parse the files for normal samples
        elif fnmatch.fnmatch(file, '*.mirna.quantification.*') and (tpat != file[13:15]):

            f = open('C:\\Users\\Sharong\\Desktop\\Level_3\\%s' % file, 'r')
            lines = f.readlines()
            #print(lines)
            for l in lines:
                line = l.split('\t')
                line= list(line)
                #if it's the file title just ignore
                if line[0].startswith("miRNA_ID"):
                    continue
                #if it's the first line of data, ignore the initial string
                elif line[0].startswith("crossed-mapped"):
                    temp = line[0][14:]
                    allcountsN[temp]= float(line[2])
                    howmanyN[temp] = 1
                    #allvalsN.append(temp)
                    #allvalsN.append(float(line[2]))
                #otherwise, just add to dictionary and keep track of sample counts
                elif line[0] in allcountsN:
                    allcountsN[line[0]]+= float(line[2])
                    howmanyN[line[0]]+=1
                    #allvalsN.append(float(line[2]))
                else:
                    allcountsN[line[0]]= float(line[2])
                    howmanyN[line[0]] = 1
                    #allvalsN.append(float(line[2]))

    #print(allcountsT)
    #print(howmanyT)
    #print("all vals:")
    #print(allvalsN)
    #print(allvalsT)

    mirna_count = 0
   #normalize:
    for k,v in allcountsT.items():
        temp = float(howmanyT[k])
        allcountsT[k]= v/temp
    for k,v in allcountsN.items():
        temp = float(howmanyN[k])
        allcountsN[k]= v/temp

    mirna_count = howmanyN[0]

    print("Tumor miRNA counts:")
    print(allcountsT)
    print("Normal miRNA counts:")
    print(allcountsN)

    return allcountsT,allcountsN,mirna_count


#sample=t if tumor, =n if normal tissue
def STDEV(dictionary,sample):
    #each miRNA key will have standard deviation value
    stdev_dict = dict()
    for k,v in dictionary.items():
        avg = dictionary[k]
        total = 0





def paired_tTest(total_mirna, allcountsN, allcountsT):
    n = total_mirna
    #store t-scores
    tscores = dict()

    #allcounts dictionaries contain the average values per miRNA
    for key,val in allcountsT.items():
        for k,v in allcountsN.items():
            continue


"""
def pval():
    #degress of freedom= n-1 = 2-1 (tumor vs. normal)
    dof = 1
"""

def main():
    tumorDict, normalDict, total_mirna = getFiles()



if __name__ == '__main__':
    main()
