

import mirnatarget
import mrnadata
import numpy as np


# for DAVID API:
import urllib2



__author__ = 'Sharong'
#import numpy as np
#from scipy.stats import ttest_1samp, wilcoxon, ttest_ind, mannwhitneyu
import csv
import glob

import os
import fnmatch
from scipy import stats
import re

#PATH = "C:\Users\Sharong\Desktop\Backed up Data 12.14.2012\Desktop\UCSD\2013-2014 academic year\BIMM 185\DB\BRCA\miRNASeq\BCGSC__IlluminaHiSeq_miRNASeq\Level_3"


def getFiles_miRNA():
    #print("Location of files: ")
    #directory = input()
    #directory = '//home//sharon//Desktop//TCGA//BRCA//miRNA'
    #directory = '//home//sharon//Desktop//TCGA//LIHC//miRNA'
    #directory = '//home//sharon//Desktop//TCGA//LUAD//miRNA'
    #directory = '//home//sharon//Desktop//TCGA//ESCA//miRNA'
    #directory = '//home//sharon//Desktop//TCGA//HNSC//miRNA'
    directory = '//home//sharon//Desktop//TCGA//KICH//miRNA'

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

            #f = open('//home//sharon//Desktop//TCGA//BRCA//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//LIHC//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//LUAD//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//ESCA//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//HNSC//miRNA//%s' % file, 'r')
            f = open('//home//sharon//Desktop//TCGA//KICH//miRNA//%s' % file, 'r')

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
                    #if RPM< .125 set to 0
                    if float(line[2]) <0.125:
                        allcountsT[line[0]]+= float(0.00000)
                        howmanyT[line[0]]+=1
                        allsamplesT[line[0]].append(float(0.00000))
                    else:
                        allcountsT[line[0]]+= float(line[2])
                        howmanyT[line[0]]+=1
                        allsamplesT[line[0]].append(float(line[2]))
                else:
                    if float(line[2]) < 0.125:
                        allcountsT[line[0]]= float(line[2])
                        howmanyT[line[0]] = 1
                        allsamplesT.setdefault(line[0],[])
                        allsamplesT[line[0]].append(float(line[2]))
                    else:
                        allcountsT[line[0]]= float(line[2])
                        howmanyT[line[0]] = 1
                        allsamplesT.setdefault(line[0],[])
                        allsamplesT[line[0]].append(float(line[2]))

        #parse the files for normal samples
        elif fnmatch.fnmatch(file, '*.mirna.quantification.txt') and (tpat != file[13:15]):


            #f = open('//home//sharon//Desktop//TCGA//BRCA//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//LIHC//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//LUAD//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//ESCA//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//HNSC//miRNA//%s' % file, 'r')
            f = open('//home//sharon//Desktop//TCGA//KICH//miRNA//%s' % file, 'r')

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
                    #if RPM < 0.125 set to 0
                    if float(line[2]) < 0.125:
                        allcountsN[line[0]]+= float(0.00000)
                        howmanyN[line[0]]+=1
                        allsamplesN[line[0]].append(float(0.00000))
                    else:
                       allcountsN[line[0]]+= float(line[2])
                       howmanyN[line[0]]+=1
                       allsamplesN[line[0]].append(float(line[2]))
                else:
                    if float(line[2]) < 0.125:
                        allcountsN[line[0]]= float(0.00000)
                        howmanyN[line[0]] = 1
                        allsamplesN.setdefault(line[0],[])
                        allsamplesN[line[0]].append(float(0.00000))
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
    print "miRNA count: ", mirna_count

    #print "tempcount tumor, normal: "
    #print tempcount_tumor, tempcount_normal

    #print "dict size: ", len(allcountsN), len(allsamplesT), len(allcountsT)

    #deleteDuplicates(allsamplesN,allsamplesT)

    #allcounts have the total RPM for that miRNA divided by the number of total samples
    #allsamples includes all the RPM values collected per each miRNA
    return allsamplesN,allsamplesT


def getfiles_mrna():

    #match between files names and barcodes
    #filemap = open('//home//sharon//Desktop//TCGA//BRCA//mRNA//file_manifest.txt','r')
    #filemap = open('//home//sharon//Desktop//TCGA//LIHC//mRNA//file_manifest.txt','r')
    #filemap = open('//home//sharon//Desktop//TCGA//LUAD//mRNA//file_manifest.txt','r')
    #filemap = open('//home//sharon//Desktop//TCGA//ESCA//mRNA//file_manifest.txt','r')
    #filemap = open('//home//sharon//Desktop//TCGA//HNSC//mRNA//file_manifest.txt','r')
    filemap = open('//home//sharon//Desktop//TCGA//KICH//mRNA//file_manifest.txt','r')


    filemap_dict = dict()
    for line in filemap:
        if line in ['\n', '\r\n']:
            continue
        else:
            line = line.strip().split()

            if line[0].startswith("Platform"):
                continue
            else:
             filemap_dict[line[6]] = line[5]

    #print filemap_dict


    #directory = '//home//sharon//Desktop//TCGA//BRCA//mRNA'
    #directory = '//home//sharon//Desktop//TCGA//LIHC//mRNA'
    #directory = '//home//sharon//Desktop//TCGA//LUAD//mRNA'
    #directory = '//home//sharon//Desktop//TCGA//ESCA//mRNA'
    #directory = '//home//sharon//Desktop//TCGA//HNSC//mRNA'
    directory = '//home//sharon//Desktop//TCGA//KICH//mRNA'


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
        #print file
        #parse files for tumor samples
        if fnmatch.fnmatch(file, '*.rsem.genes.normalized_results') and (tpat==filemap_dict[file][13:15]):
            #print(file[13:15])
            #f = open('//home//sharon//Desktop//TCGA//BRCA//mRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//LIHC//mRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//LUAD//mRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//ESCA//mRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//HNSC//mRNA//%s' % file, 'r')
            f = open('//home//sharon//Desktop//TCGA//KICH//mRNA//%s' % file, 'r')


            lines = f.readlines()

            for l in lines:
                line = l.split('\t')

                #if it's the file title just ignore
                if line[0].startswith("gene_id") or line[0].startswith("?"):
                    continue

                else:
                    #id[0] will be gene name, id[1] will be gene id #
                    id = line[0].split('|')
                    #id[1]= id[1][1:]
                    norm_count = line[1]
                    id[0] = id[0].lower()



                    #otherwise, just add to dictionary and keep track of sample counts
                    if id[0] in allcountsT:
                        #if norm_count< .125 set to 0
                        if float(norm_count) <0.125:
                            allcountsT[id[0]]+= float(0.00000)
                            howmanyT[id[0]]+=1
                            allsamplesT[id[0]].append(float(0.00000))
                        else:
                            allcountsT[id[0]]+= float(norm_count)
                            howmanyT[id[0]]+=1
                            allsamplesT[id[0]].append(float(norm_count))
                    else:
                        if float(norm_count) < 0.125:
                            allcountsT[id[0]]= float(norm_count)
                            howmanyT[id[0]] = 1
                            allsamplesT.setdefault(id[0],[])
                            allsamplesT[id[0]].append(float(norm_count))
                        else:
                            allcountsT[id[0]]= float(norm_count)
                            howmanyT[id[0]] = 1
                            allsamplesT.setdefault(id[0],[])
                            allsamplesT[id[0]].append(float(norm_count))
        elif fnmatch.fnmatch(file, '*.rsem.genes.normalized_results') and ('11'==filemap_dict[file][13:15]):
            #print(file[13:15])
            #f = open('//home//sharon//Desktop//TCGA//BRCA//mRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//LIHC//mRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//LUAD//mRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//ESCA//mRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//HNSC//mRNA//%s' % file, 'r')
            f = open('//home//sharon//Desktop//TCGA//KICH//mRNA//%s' % file, 'r')

            lines = f.readlines()

            for l in lines:
                line = l.split('\t')

                #if it's the file title just ignore
                if line[0].startswith("gene_id") or line[0].startswith("?"):
                    continue

                else:
                    #id[0] will be gene name, id[1] will be gene id #
                    id = line[0].split('|')
                    #id[1]= id[1][1:]
                    norm_count = line[1]
                    id[0] = id[0].lower()

                    #otherwise, just add to dictionary and keep track of sample counts
                    if id[0] in allcountsN:
                        #if norm_count< .125 set to 0
                        if float(norm_count) <0.125:
                            allcountsN[id[0]]+= float(0.00000)
                            howmanyN[id[0]]+=1
                            allsamplesN[id[0]].append(float(0.00000))
                        else:
                            allcountsN[id[0]]+= float(norm_count)
                            howmanyN[id[0]]+=1
                            allsamplesN[id[0]].append(float(norm_count))
                    else:
                        if float(norm_count) < 0.125:
                            allcountsN[id[0]]= float(0.00000)
                            howmanyN[id[0]] = 1
                            allsamplesN.setdefault(id[0],[])
                            allsamplesN[id[0]].append(float(0.00000))
                        else:
                            allcountsN[id[0]]= float(norm_count)
                            howmanyN[id[0]] = 1
                            allsamplesN.setdefault(id[0],[])
                            allsamplesN[id[0]].append(float(norm_count))

    gene_count = len(allsamplesT)
    #print gene_count


    return allsamplesT, allsamplesN, gene_count


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




def foldChange(allsampleN, allsamplesT):
    foldchange_dict = dict()
    result = dict()
    for (k,v), (k2,v2) in zip(allsamplesT.items(),allsampleN.items()):
        #normal samples
        n = 0
        for i in v:
            n+= i
        #tumor samples
        t = 0
        for j in v2:
            t+= j
        if n==0:
            foldchange_dict[k] = float(0)
        else:
            foldchange_dict[k] = np.log2(float((t+0.1) / (n+0.1)))
        if foldchange_dict[k] >=1.0 or foldchange_dict[k]<=-1.0:
            result[k] = foldchange_dict[k]
            #print foldchange_dict[k]
    #print "fold change:"
    #print(foldchange_dict)
    return result




def paired_tTest(total_mirna, allsamplesN, allsamplesT):

    #output of all differentially expressed miRNAs, include amounts?
    #output = open('//home//sharon//Desktop//TCGA//BRCA//output_miRNA.txt', 'w')
    #output = open('//home//sharon//Desktop//TCGA//LIHC//output_miRNA.txt', 'w')
    #output = open('//home//sharon//Desktop//TCGA//LUAD//output_miRNA.txt', 'w')
    #output = open('//home//sharon//Desktop//TCGA//ESCA//output_miRNA.txt', 'w')
    #output = open('//home//sharon//Desktop//TCGA//HNSC//output_miRNA.txt', 'w')
    output = open('//home//sharon//Desktop//TCGA//KICH//output_miRNA.txt', 'w')

    output.writelines("miRNA\tlog2_FC\n")

    #stores all the significant p-vals <0.05
    pval_dict = dict()

    FC = foldChange(allsamplesN,allsamplesT)

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

            #signficant if p<0.05 and FC = +-1
            if paired_sample[1] < 0.05 and k in FC:
                if (FC[k]>=1 or FC[k]<=-1 ):
                    #print("paired sample:", paired_sample)
                    significant_count +=1
                    pval_dict[k] = paired_sample[1]


                    output.write(k + "\t"+ str(FC[k])+"\n")
    print "significant count: ", significant_count
    output.close()
    return pval_dict


#prepares the layout for XGMML output file for Cytoscape
# Gene_ID, gene_name, and then each column has expression and p-values columns
def outputForCyto(targets, pval_dict):
    #add log2 FC and P-value to dictionary as values
    #miRNA_FC = open('//home//sharon//Desktop//TCGA//BRCA//output.txt', 'r')
    #miRNA_FC = open('//home//sharon//Desktop//TCGA//LIHC//output.txt', 'r')
    #miRNA_FC = open('//home//sharon//Desktop//TCGA//LUAD//output.txt', 'r')
    #miRNA_FC = open('//home//sharon//Desktop//TCGA//ESCA//output.txt', 'r')
    #miRNA_FC = open('//home//sharon//Desktop//TCGA//HNSC//output.txt', 'r')
    miRNA_FC = open('//home//sharon//Desktop//TCGA//KICH//output.txt', 'r')

    #cyto_output = open('//home//sharon//Desktop//TCGA//BRCA//output_cyto.txt', 'w+')
    #cyto_output = open('//home//sharon//Desktop//TCGA//LIHC//output_cyto.txt', 'w+')
    #cyto_output = open('//home//sharon//Desktop//TCGA//LUAD//output_cyto.txt', 'w+')
    #cyto_output = open('//home//sharon//Desktop//TCGA//ESCA//output_cyto.txt', 'w+')
    #cyto_output = open('//home//sharon//Desktop//TCGA//HNSC//output_cyto.txt', 'w+')
    cyto_output = open('//home//sharon//Desktop//TCGA//KICH//output_cyto.txt', 'w+')

    #genelist_output = open('//home//sharon//Desktop//TCGA//BRCA//output_genelist.txt', 'w+')
    #genelist_output = open('//home//sharon//Desktop//TCGA//LIHC//output_genelist.txt', 'w+')
    #genelist_output = open('//home//sharon//Desktop//TCGA//LUAD//output_genelist.txt', 'w+')
    #genelist_output = open('//home//sharon//Desktop//TCGA//ESCA//output_genelist.txt', 'w+')
    #genelist_output = open('//home//sharon//Desktop//TCGA//HNSC//output_genelist.txt', 'w+')
    genelist_output = open('//home//sharon//Desktop//TCGA//KICH//output_genelist.txt', 'w+')

    cyto_output.write("gene_id\tgene_name\tlog2_fold_change\tp_value\n")
    miRNA_FC_dict = dict()

    #keep track of all genes already written to file
    genes = []

    for l in miRNA_FC:
        l = l.strip().split()
        miRNA_FC_dict[l[0]] = float(l[1])
    #print miRNA_FC_dict
    #print targets
    for k,v in targets.items():
        #print k[0] , miRNA_FC_dict[k]
        if k[0] in miRNA_FC_dict:
            if str(k[2]) not in genes:
                genes.append(str(k[2]))
                genelist_output.write(k[2].lower()+"\n")
            cyto_output.write(k[1]+"\t"+str(k[2])+"\t"+str(miRNA_FC_dict[k[0]])+"\t"+str(pval_dict[k[0]])+"\n")
    cyto_output.close()
    miRNA_FC.close()


def outputForCyto2(targets, pval_dict):
    #add log2 FC and P-value to dictionary as values

    #miRNA_FC = open('//home//sharon//Desktop//TCGA//BRCA//output_miRNA.txt', 'r')
    #miRNA_FC = open('//home//sharon//Desktop//TCGA//LIHC//output_miRNA.txt', 'r')
    #miRNA_FC = open('//home//sharon//Desktop//TCGA//LUAD//output_miRNA.txt', 'r')
    #miRNA_FC = open('//home//sharon//Desktop//TCGA//ESCA//output_miRNA.txt', 'r')
    #miRNA_FC = open('//home//sharon//Desktop//TCGA//HNSC//output_miRNA.txt', 'r')
    miRNA_FC = open('//home//sharon//Desktop//TCGA//KICH//output_miRNA.txt', 'r')

    #cyto_output = open('//home//sharon//Desktop//TCGA//BRCA//output_cyto.txt', 'w')
    #cyto_output = open('//home//sharon//Desktop//TCGA//LIHC//output_cyto.txt', 'w')
    #cyto_output = open('//home//sharon//Desktop//TCGA//LUAD//output_cyto.txt', 'w')
    #cyto_output = open('//home//sharon//Desktop//TCGA//ESCA//output_cyto.txt', 'w')
    #cyto_output = open('//home//sharon//Desktop//TCGA//HNSC//output_cyto.txt', 'w')
    cyto_output = open('//home//sharon//Desktop//TCGA//KICH//output_cyto.txt', 'w')

    #genelist_output = open('//home//sharon//Desktop//TCGA//BRCA//output_genelist.txt', 'w')
    #genelist_output = open('//home//sharon//Desktop//TCGA//LIHC//output_genelist.txt', 'w')
    #genelist_output = open('//home//sharon//Desktop//TCGA//LUAD//output_genelist.txt', 'w')
    #genelist_output = open('//home//sharon//Desktop//TCGA//ESCA//output_genelist.txt', 'w')
    #genelist_output = open('//home//sharon//Desktop//TCGA//HNSC//output_genelist.txt', 'w')
    genelist_output = open('//home//sharon//Desktop//TCGA//KICH//output_genelist.txt', 'w')

    cyto_output.writelines("#differentially expressed miRNA target genes. FC threshold -+1, p<0.05"+ "\n")
    cyto_output.writelines("gene_name\tlog2_fold_change\tp_value\n")

    genelist_output.writelines("#miRNA targets of differentially expressed miRNAs from both CLASH and mirTarBase\n")

    miRNA_FC_dict = dict()

    #keep track of all genes already written to file
    genes = []

    for l in miRNA_FC:
        if l.startswith("miRNA"):
            continue
        else:
            l = l.strip().split()
            miRNA_FC_dict[l[0]] = float(l[1])
    #print miRNA_FC_dict
    #print targets
    for k,v in targets.items():
        if k in miRNA_FC_dict:
            if str(v) not in genes:
                genes.append(str(v))
                genelist_output.write(v.lower()+"\n")
            cyto_output.write(v.lower()+"\t"+str(miRNA_FC_dict[k])+"\t"+str(pval_dict[k])+"\n")
    cyto_output.close()
    miRNA_FC.close()



def outputForCyto3(targets, pval_dict):
    #add log2 FC and P-value to dictionary as values

    #miRNA_FC = open('//home//sharon//Desktop//TCGA//BRCA//output_miRNA.txt', 'r')
    #miRNA_FC = open('//home//sharon//Desktop//TCGA//LIHC//output_miRNA.txt', 'r')
    #miRNA_FC = open('//home//sharon//Desktop//TCGA//LUAD//output_miRNA.txt', 'r')
    #miRNA_FC = open('//home//sharon//Desktop//TCGA//ESCA//output_miRNA.txt', 'r')
    #miRNA_FC = open('//home//sharon//Desktop//TCGA//HNSC//output_miRNA.txt', 'r')
    miRNA_FC = open('//home//sharon//Desktop//TCGA//KICH//output_miRNA.txt', 'r')

    #cyto_output = open('//home//sharon//Desktop//TCGA//BRCA//output_cyto.txt', 'a')
    #cyto_output = open('//home//sharon//Desktop//TCGA//LIHC//output_cyto.txt', 'a')
    #cyto_output = open('//home//sharon//Desktop//TCGA//LUAD//output_cyto.txt', 'a')
    #cyto_output = open('//home//sharon//Desktop//TCGA//ESCA//output_cyto.txt', 'a')
    #cyto_output = open('//home//sharon//Desktop//TCGA//HNSC//output_cyto.txt', 'a')
    cyto_output = open('//home//sharon//Desktop//TCGA//KICH//output_cyto.txt', 'a')

    #genelist_output = open('//home//sharon//Desktop//TCGA//BRCA//output_genelist.txt', 'a+')
    #genelist_output = open('//home//sharon//Desktop//TCGA//LIHC//output_genelist.txt', 'a+')
    #genelist_output = open('//home//sharon//Desktop//TCGA//LUAD//output_genelist.txt', 'a+')
    #genelist_output = open('//home//sharon//Desktop//TCGA//ESCA//output_genelist.txt', 'a+')
    #genelist_output = open('//home//sharon//Desktop//TCGA//HNSC//output_genelist.txt', 'a+')
    genelist_output = open('//home//sharon//Desktop//TCGA//KICH//output_genelist.txt', 'a+')

    #cyto_output.write("gene_name\tmiRNA_target_log2_FC\tp_value\n")
    miRNA_FC_dict = dict()

    #keep track of all genes already written to file
    genes = []

    #since mirTarBase is called after CLASH database, first see which genes already found
    for line in genelist_output:
        line = line.strip().lower()
        if line.startswith("#"):
            continue
        else:
            genes.append(line)
    print "genelist: "
    print genes

    for l in miRNA_FC:
        if l.startswith("miRNA"):
            continue
        else:
            l = l.strip().split()
            miRNA_FC_dict[l[0]] = float(l[1])
    print miRNA_FC_dict
    print targets
    for k,v in targets.items():
        #print v[0]
        if k in miRNA_FC_dict:
            if str(v[0]) not in genes:
                genes.append(str(v[0]))
                genelist_output.write(v[0].lower()+"\n")
            cyto_output.write(v[0].lower()+"\t"+str(miRNA_FC_dict[k])+"\t"+str(pval_dict[k])+"\n")
    cyto_output.close()
    miRNA_FC.close()

"""
API:
http://david.abcc.ncifcrf.gov/api.jsp?type=xxxxx&ids=XXXXX,XXXXX,XXXXXX,&tool=xxxx&annot=xxxxx,xxxxxx,xxxxx,
    type  =  one of DAVID recognized gene types
    annot  = a list of desired annotation  categories separated by ","
    ids  = a list of user's gene IDs separated by ","
    tool  = one of DAVID tool names

"""
def DAVID_annotation():
    IDs = ""
    #genelist_output = open('//home//sharon//Desktop//TCGA//BRCA//output_genelist.txt', 'r')
    #genelist_output = open('//home//sharon//Desktop//TCGA//LIHC//output_genelist.txt', 'r')
    #genelist_output = open('//home//sharon//Desktop//TCGA//LUAD//output_genelist.txt', 'r')
    #genelist_output = open('//home//sharon//Desktop//TCGA//ESCA//output_genelist.txt', 'r')
    #genelist_output = open('//home//sharon//Desktop//TCGA//HNSC//output_genelist.txt', 'r')
    genelist_output = open('//home//sharon//Desktop//TCGA//KICH//output_genelist.txt', 'r')


    for line in genelist_output:
        line = line.strip()
        IDs+=line+","
    print IDs
    url = "http://david.abcc.ncifcrf.gov/api.jsp?type=GENE_SYMBOL&ids="+IDs+"&tool=term2term&annot=GOTERM_BP_FAT"
    try:
      result = urllib2.urlopen(url)
      print result.read()
    except urllib2.URLError, e:
        print e.reason


#check for intersection between gene lists from miRNA and mRNA
def combine_data():
    #miRNA_genes = open('//home//sharon//Desktop//TCGA//BRCA//output_genelist.txt','r')
    #miRNA_genes = open('//home//sharon//Desktop//TCGA//LIHC//output_genelist.txt','r')
    #miRNA_genes = open('//home//sharon//Desktop//TCGA//LUAD//output_genelist.txt','r')
    #miRNA_genes = open('//home//sharon//Desktop//TCGA//ESCA//output_genelist.txt','r')
    #miRNA_genes = open('//home//sharon//Desktop//TCGA//HNSC//output_genelist.txt','r')
    miRNA_genes = open('//home//sharon//Desktop//TCGA//KICH//output_genelist.txt','r')

    #mRNA_gnes = open("//home//sharon//Desktop//TCGA//BRCA//output_mRNA_genes.txt",'r')
    #mRNA_gnes = open("//home//sharon//Desktop//TCGA//LIHC//output_mRNA_genes.txt",'r')
    #mRNA_gnes = open("//home//sharon//Desktop//TCGA//LUAD//output_mRNA_genes.txt",'r')
    #mRNA_gnes = open("//home//sharon//Desktop//TCGA//ESCA//output_mRNA_genes.txt",'r')
    #mRNA_gnes = open("//home//sharon//Desktop//TCGA//HNSC//output_mRNA_genes.txt",'r')
    mRNA_gnes = open("//home//sharon//Desktop//TCGA//KICH//output_mRNA_genes.txt",'r')

    miRNA_list = []
    mRNA_list = []
    mRNA_dict = dict()

    for line in miRNA_genes:
        line = line.strip()
        if line.startswith("gene_name") or line.startswith("#"):
            continue
        else:
            miRNA_list.append(line)
    for line in mRNA_gnes:
        if line.startswith("gene_name") or line.startswith("#"):
            continue
        elif line == '\n':
            continue
        else:
            line = line.split()
            mRNA_list.append(line[0])

            #if no gene_ID available:
            if len(line) == 2:
                mRNA_dict[(line[0],'N/A')] = line[1]
            else:
                #key: gene_name, gene_ID ; value: FC
                mRNA_dict[(line[0],line[1])] = line[2]

    #biomarkers = open('//home//sharon//Desktop//TCGA//BRCA//output_biomarkers.txt','w')
    #biomarkers = open('//home//sharon//Desktop//TCGA//LIHC//output_biomarkers.txt','w')
    #biomarkers = open('//home//sharon//Desktop//TCGA//LUAD//output_biomarkers.txt','w')
    #biomarkers = open('//home//sharon//Desktop//TCGA//ESCA//output_biomarkers.txt','w')
    #biomarkers = open('//home//sharon//Desktop//TCGA//HNSC//output_biomarkers.txt','w')
    biomarkers = open('//home//sharon//Desktop//TCGA//KICH//output_biomarkers.txt','w')


    biomarkers.write("#Genes differentially expressed in both miRNA and mRNA data\n")

    """
    print "both lists:"
    print miRNA_list
    print mRNA_list

    for g in miRNA_list:

        if g in mRNA_list:
            print "g " + g
            biomarkers.write(g +"\n")
    """

    #output with gene from miRNA and FC from both miRNA and mRNA data sets
    #FC_both = open('//home//sharon//Desktop//TCGA//BRCA//output_combineFC.txt','w')
    #FC_both = open('//home//sharon//Desktop//TCGA//LIHC//output_combineFC.txt','w')
    #FC_both = open('//home//sharon//Desktop//TCGA//LUAD//output_combineFC.txt','w')
    #FC_both = open('//home//sharon//Desktop//TCGA//ESCA//output_combineFC.txt','w')
    #FC_both = open('//home//sharon//Desktop//TCGA//HNSC//output_combineFC.txt','w')
    FC_both = open('//home//sharon//Desktop//TCGA//KICH//output_combineFC.txt','w')

    FC_both.write("#Genes differentially expressed in both miRNA and mRNA data\n")
    FC_both.write("gene_name\tgene_ID\tlog2_FC_miRNA\tlog2_FC_mRNA\n")


    #cyto_output = open('//home//sharon//Desktop//TCGA//BRCA//output_cyto.txt', 'r')
    #cyto_output = open('//home//sharon//Desktop//TCGA//LIHC//output_cyto.txt', 'r')
    #cyto_output = open('//home//sharon//Desktop//TCGA//LUAD//output_cyto.txt', 'r')
    #cyto_output = open('//home//sharon//Desktop//TCGA//ESCA//output_cyto.txt', 'r')
    #cyto_output = open('//home//sharon//Desktop//TCGA//HNSC//output_cyto.txt', 'r')
    cyto_output = open('//home//sharon//Desktop//TCGA//KICH//output_cyto.txt', 'r')

    cyto = dict()
    for line in cyto_output:

        line = line.strip()
        if line.startswith("gene_name") or line.startswith("#"):
            continue
        else:
            l = re.split(r'\s',line)
            # key= gene_name ; values= FC, p_val
            cyto[l[0]] = (l[1],l[2])


    for (k,v),(k2,v2) in zip(cyto.items(),mRNA_dict.items()):
        FC_both.write(k2[0]+"\t"+k2[1]+"\t"+v[0]+"\t"+v2+"\n")
        biomarkers.writelines(str(k)+"\n")

    FC_both.close()
    cyto_output.close()



#Pearson and Spearman
#S depicts monotonic relationships while P depicts linear relationships
def correlation():
    #output with gene from miRNA and FC from both miRNA and mRNA data sets
    #FC_both = open('//home//sharon//Desktop//TCGA//BRCA//output_combineFC.txt','r')
    #FC_both = open('//home//sharon//Desktop//TCGA//LIHC//output_combineFC.txt','r')
    #FC_both = open('//home//sharon//Desktop//TCGA//LUAD//output_combineFC.txt','r')
    #FC_both = open('//home//sharon//Desktop//TCGA//ESCA//output_combineFC.txt','r')
    #FC_both = open('//home//sharon//Desktop//TCGA//HNSC//output_combineFC.txt','r')
    FC_both = open('//home//sharon//Desktop//TCGA//KICH//output_combineFC.txt','r')

    genename = []
    #mirna FC:
    x = []
    #mrna FC:
    y = []

    for line in FC_both:
        line = line.strip()
        if line.startswith("#Genes") or line.startswith("gene_name"):
            continue
        else:
            line = line.split()
            genename.append(line[0])
            if "inf" in line[2]:
                if line[2].startswith("-"):
                    x.append(float(-100000.00))
                else:
                    x.append(float(100000.00))
            else:
                x.append(float(line[2]))
            if "inf" in line[3]:
                if line[3].startswith("-"):
                    y.append(float(-100000.00))
                else:
                    y.append(float(100000.00))
            else:
                y.append(float(line[3]))


            #print x
            #print y

    #pearson:
    result_pearson =  stats.pearsonr(x,y)
    print(result_pearson)

    #spearman:
    result_spearman =  stats.spearmanr(x,y)
    print(result_spearman)



def miRNAtargets(targets_clash,targets_mirtar,allsamplesN_mirna,allsamplesT_mirna):

    #process all miRNA targets to once dictionary
    genetargets = dict()

    for k,v in targets_clash.items():
        if k in genetargets:
            if v in genetargets[k]:
                continue
            else:
                genetargets[k].append(v)
        else:
            genetargets.setdefault(k,[])
            genetargets[k].append(v)
    for k,v in targets_mirtar.items():
        if k in genetargets:
            if v[0] in genetargets[k]:
                continue
            else:
                genetargets[k].append(v[0])
        else:
            genetargets.setdefault(k,[])
            genetargets[k].append(v[0])

    #print "genetargets:"
    #print genetargets

     #tumor
    alltargetsT = dict()
    for k,v in allsamplesT_mirna.items():
        mirna = str(k)

        v = str(v).replace("[","").replace("]","")

        if mirna in genetargets:
            if mirna in alltargetsT:
                continue
            else:
                l=1
                tmp = str(genetargets[mirna]).replace("[","").replace("]","")
                tmp = tmp.replace("'","")
                if "," in tmp:
                    tmp = tmp.split(",")
                    l = len(tmp)
                #tmp = tmp.split(",")
                #tmp = str(genetargets[mirna])[2:-2]
                #print tmp
                #l = len(genetargets[mirna])



                #print k, v
                if l==1:
                    #print "l=1",patient,tmp,v
                    alltargetsT[(tmp).strip()] = v
                elif l>1:
                    for i in range(l):
                        alltargetsT[(tmp[i]).strip()] = v
                        #print "yay ",patient, tmp

    #print alltargetsT


    #normal
    alltargetsN= dict()
    for k,v in allsamplesN_mirna.items():
        mirna = str(k)

        v = str(v).replace("[","").replace("]","")

        if mirna in genetargets:
            l=1
            tmp = str(genetargets[mirna]).replace("[","").replace("]","")
            tmp = tmp.replace("'","")
            if "," in tmp:
                tmp = tmp.split(",")
                l = len(tmp)
            #tmp = tmp.split(",")
            #tmp = str(genetargets[mirna])[2:-2]
            #print tmp
            #l = len(genetargets[mirna])

            #print l
            if l==1:
                #print "l=1",patient,tmp,v
                alltargetsN[(tmp).strip()] = v
            elif l>1:
                for i in range(l):
                    alltargetsN[(tmp[i]).strip()] = v
                    #print "yay ",patient, tmp
    #print "normal targets:"
    #print alltargetsN


    FC_each = dict()
    for (k,v),(k2,v2) in zip(alltargetsT.items(),alltargetsN.items()):
        if k==k2:
            tmp = []
            v = [float(x) for x in v.split(",")]
            v2 = [float(x) for x in v2.split(",")]
            for i in range(0,len(v)):
                #print len(v), v
                tmp.append((float(v[i])+0.1)/((float(v2[i]))+0.1))
            FC_each[k] = tmp


    print FC_each
    return FC_each

def pair_mrna(allsamplesT_mrna,allsamplesN_mrna):
    FC_each = dict()
    for (k,v),(k2,v2) in zip(allsamplesT_mrna.items(),allsamplesN_mrna.items()):
        if k==k2:
            tmp = []
            v = [float(x) for x in v]
            v2 = [float(x) for x in v2]
            for i in range(0,len(v)):
                #print len(v), v
                tmp.append((float(v[i])+0.1)/((float(v2[i]))+0.1))
            FC_each[k] = tmp
    print FC_each
    return FC_each


def mutual_genes(FC_each_mrna,FC_each_mirna):

    for k in FC_each_mrna.keys():
        for k2 in FC_each_mirna.keys():
            if k==k2:
                print k

def main():

    allsamplesN, allsamplesT = getFiles_miRNA()
    #pval_dict = paired_tTest(allsamplesN,allsamplesT)

    print "dicts:"
    print allsamplesT
    print allsamplesN


    #returns a dictionary: 3 keys[miRNA_ID,gene_ID,gene_name] and miSVR_score as value
    ##TEMP!! need to work with CLASH DB to narrow down predictions
    #work with DAVID annotation tool and import to Bader lab enchrichemnt map app
    #targets = mirnatarget.microRNA()
    targets_clash = mirnatarget.CLASH()
    targets_mirtar = mirnatarget.MirTarBase()

    print "targets:"
    print targets_clash
    print targets_mirtar

    FC_each_mirna = miRNAtargets(targets_clash,targets_mirtar,allsamplesN,allsamplesT)

    allsamplesT_mrna, allsamplesN_mrna, gene_count = getfiles_mrna()

    #print "mrna data:"
    #print allsamplesT_mrna
    #print allsamplesN_mrna

    FC_each_mrna = pair_mrna(allsamplesT_mrna,allsamplesN_mrna)

    mutual_genes(FC_each_mrna,FC_each_mirna)
    """
    #outputForCyto(targets, pval_dict)
    outputForCyto2(targets_clash, pval_dict)

    #add gene results from outputForCyto2:
    outputForCyto3(targets_mirtar,pval_dict)

    #DAVID_annotation()

    mrnadata.main()

    combine_data()

    correlation()
    """

if __name__ == '__main__':
    main()