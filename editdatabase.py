# TCGA file parser


__author__ = 'Sharong'
#import numpy as np
#from scipy.stats import ttest_1samp, wilcoxon, ttest_ind, mannwhitneyu
import csv
import glob

import os
import fnmatch
from scipy import stats

#PATH = "C:\Users\Sharong\Desktop\Backed up Data 12.14.2012\Desktop\UCSD\2013-2014 academic year\BIMM 185\DB\BRCA\miRNASeq\BCGSC__IlluminaHiSeq_miRNASeq\Level_3"


def getFiles():
    #print("Location of files: ")
    #directory = input()
    directory = '//home//sharon//Desktop//TCGA//BRCA//BRCA_allinone//edit'

    path = r"%s" % directory
    #dictionary for normal and tumor samples
    allcountsN = dict()
    allcountsT = dict()

    #counts how many miRNA data I have per each
    howmanyN = dict()
    howmanyT = dict()

    #01 at a specific location means tumor sample
    tpat = "01"

    #dictionary for normal and tumor samples-- keeps all values
    allsamplesN = dict()
    allsamplesT = dict()

    for file in os.listdir(path):
        #parse files for tumor samples
        if fnmatch.fnmatch(file, '*.mirna.quantification.txt') and (tpat==file[13:15]):
            #print(file[13:15])
            f = open('//home//sharon//Desktop//TCGA//BRCA//BRCA_allinone//edit//%s' % file, 'r')
            lines = f.readlines()
            #print(lines)
            for l in lines:

                line = l.split('\t')
                line= list(line)
                #if it's the file title just ignore
                if line[0].startswith("miRNA_ID"):
                    continue
                #otherwise, just add to dictionary and keep track of sample counts
                elif line[0] in allcountsT:
                    allcountsT[line[0]]+= float(line[2])
                    howmanyT[line[0]]+=1
                    allsamplesT[line[0]].append(float(line[2]))
                else:
                    allcountsT[line[0]]= float(line[2])
                    howmanyT[line[0]] = 1

                    allsamplesT.setdefault(line[0],[])
                    allsamplesT[line[0]].append(float(line[2]))

        #parse the files for normal samples
        elif fnmatch.fnmatch(file, '*.mirna.quantification.txt') and (tpat != file[13:15]):

            #f = open('C:\\Users\\Sharong\\Desktop\\Level_3\\%s' % file, 'r')
            f = open('//home//sharon//Desktop//TCGA//BRCA//BRCA_allinone//edit//%s' % file, 'r')

            lines = f.readlines()
            #print(lines)
            for l in lines:
                line = l.split('\t')
                line= list(line)
                #if it's the file title just ignore
                if line[0].startswith("miRNA_ID"):
                    continue
                #otherwise, just add to dictionary and keep track of sample counts
                elif line[0] in allcountsN:
                    allcountsN[line[0]]+= float(line[2])
                    howmanyN[line[0]]+=1
                    allsamplesN[line[0]].append(float(line[2]))
                else:
                    allcountsN[line[0]]= float(line[2])
                    howmanyN[line[0]] = 1
                    allsamplesN.setdefault(line[0],[])
                    allsamplesN[line[0]].append(float(line[2]))

    #testing:

    #for k,v in allsamplesN.items():
        #if k=='hsa-mir-1322':
            #s = sum([float(x) for x in v])
            #print("sum", s)
            #sum+= float(v)


    mirna_count = 0
   #normalize: TEMP
    for k,v in allcountsT.items():
        temp = float(howmanyT[k])
        allcountsT[k]= v/temp
    for k,v in allcountsN.items():
        temp = float(howmanyN[k])
        allcountsN[k]= v/temp


    mirna_count = len(allcountsT)

    #print "tempcount tumor, normal: "
    #print tempcount_tumor, tempcount_normal

    #print "dict size: ", len(allcountsN), len(allsamplesT), len(allcountsT)

    #deleteDuplicates(allsamplesN,allsamplesT)

    return allcountsT,allcountsN,mirna_count, allsamplesN,allsamplesT

"""
#to make sure all samples have values for all miRNAs. to do  paired t-test
def deleteDuplicates(allsamplesN, allsamplesT):
    for k,v in allsamplesN.items():
        if k not in allsamplesT:
            print "DELETED"
            del allsamplesN[k]
        elif len(v) > len(allsamplesT[k]):
            print "NOT EQUAL LENGTHS!!"
    for k,v in allsamplesT.items():
        if k not in allsamplesN:
            print "DELETED"
            del allsamplesT[k]
"""


# z = (expression in tumor sample) - (mean expression in normal sample)/ (standard deviation of expression in normal sample)
def zScore(allsamplesN, allsamplesT):
    import numpy
    mean_dict = dict()
    stdev_dict = dict()
    for k,v in allsamplesN.items():
        mean_dict[k] = numpy.mean(v)
        stdev_dict[k] = numpy.std(v)

    z_score = dict()
    for k,v in allsamplesT.items():
        z_score.setdefault(k,[])
        for i in v:
            z_score[k].append((i - mean_dict[k]) / stdev_dict[k])
    print z_score



def paired_tTest(total_mirna, allsamplesN, allsamplesT):

    n = total_mirna
    significant_count = 0
    #checked - all miRNA collected in all samples
    #print("lengths: ", allsamplesT.__sizeof__(), allsamplesN.__sizeof__())
    for k,v in allsamplesT.items():
        #number of miRNAs apparently don't match so only use the ones that do
        if k in allsamplesN:

            normal = allsamplesN[k]
            tumor = allsamplesT[k]

            paired_sample = stats.ttest_rel(normal,tumor)


            if paired_sample[1] < 0.05:
                #print("paired sample:", paired_sample)
                significant_count +=1
    print "significant count: ", significant_count






def main():
    tumorDict, normalDict, total_mirna, allsamplesN, allsamplesT = getFiles()
    paired_tTest(total_mirna,allsamplesN,allsamplesT)
    zScore(allsamplesN,allsamplesT)



if __name__ == '__main__':
    main()
