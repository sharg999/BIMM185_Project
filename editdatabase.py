# TCGA file parser

import mirnatarget
import mrnadata
import numpy as np



__author__ = 'Sharong'

import os
import fnmatch
from scipy import stats
import re

def getFiles_miRNA():

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

            #f = open('//home//sharon//Desktop//TCGA//BRCA//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//LIHC//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//LUAD//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//ESCA//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//HNSC//miRNA//%s' % file, 'r')
            f = open('//home//sharon//Desktop//TCGA//KICH//miRNA//%s' % file, 'r')

            lines = f.readlines()
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
                        allcountsT[line[0]]= float(0.0000)
                        howmanyT[line[0]] = 1
                        allsamplesT.setdefault(line[0],[])
                        allsamplesT[line[0]].append(float(0.0000))
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
                        allcountsN[line[0]]= float(0.0000)
                        howmanyN[line[0]] = 1
                        allsamplesN.setdefault(line[0],[])
                        allsamplesN[line[0]].append(float(0.0000))
                    else:
                        allcountsN[line[0]]= float(line[2])
                        howmanyN[line[0]] = 1
                        allsamplesN.setdefault(line[0],[])
                        allsamplesN[line[0]].append(float(line[2]))


    mirna_count = 0

    for k,v in allcountsT.items():
        temp = float(howmanyT[k])
        allcountsT[k]= v/temp
    for k,v in allcountsN.items():
        temp = float(howmanyN[k])
        allcountsN[k]= v/temp


    mirna_count = len(allcountsT)
    print "Total miRNAs analyzed: ", mirna_count

    #allcounts have the total RPM for that miRNA divided by the number of total samples
    #allsamples includes all the RPM values collected per each miRNA
    return allcountsT,allcountsN,mirna_count, allsamplesN,allsamplesT


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
    return result




def paired_tTest(total_mirna, allsamplesN, allsamplesT):

    #output of all differentially expressed miRNAs
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
    for k,v in allsamplesT.items():
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
    print "Number of significant differentially expressed miRNAs: ", significant_count
    output.close()
    return pval_dict



#differentially expressed miRNA target genes. FC threshold -+1, p<0.05
#gene_name	log2_fold_change	p_value
def outputFCmirna1(targets, pval_dict):

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
    print "writing differentially expressed miRNA target results in output_cyto.txt..."

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
    for k,v in targets.items():
        if k in miRNA_FC_dict:
            if str(v) not in genes:
                genes.append(str(v))
                genelist_output.write(v.lower()+"\n")
            cyto_output.write(v.lower()+"\t"+str(miRNA_FC_dict[k])+"\t"+str(pval_dict[k])+"\n")
    cyto_output.close()
    miRNA_FC.close()



def outputFCmirna2(targets, pval_dict):

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

    for l in miRNA_FC:
        if l.startswith("miRNA"):
            continue
        else:
            l = l.strip().split()
            miRNA_FC_dict[l[0]] = float(l[1])

    for k,v in targets.items():
        if k in miRNA_FC_dict:
            if str(v[0]) not in genes:
                genes.append(str(v[0]))
                genelist_output.write(v[0].lower()+"\n")
            cyto_output.write(v[0].lower()+"\t"+str(miRNA_FC_dict[k])+"\t"+str(pval_dict[k])+"\n")
    cyto_output.close()
    miRNA_FC.close()


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
    print "writing mutual mRNA-miRNA differentially expressed gene results in output_biomarkers.txt..."


    biomarkers.write("#Genes differentially expressed in both miRNA and mRNA data\n")

    #output with gene from miRNA and FC from both miRNA and mRNA data sets
    #FC_both = open('//home//sharon//Desktop//TCGA//BRCA//output_combineFC.txt','w')
    #FC_both = open('//home//sharon//Desktop//TCGA//LIHC//output_combineFC.txt','w')
    #FC_both = open('//home//sharon//Desktop//TCGA//LUAD//output_combineFC.txt','w')
    #FC_both = open('//home//sharon//Desktop//TCGA//ESCA//output_combineFC.txt','w')
    #FC_both = open('//home//sharon//Desktop//TCGA//HNSC//output_combineFC.txt','w')
    FC_both = open('//home//sharon//Desktop//TCGA//KICH//output_combineFC.txt','w')

    FC_both.write("#Genes differentially expressed in both miRNA and mRNA data\n")
    FC_both.write("gene_name\tgene_ID\tlog2_FC_miRNA\tlog2_FC_mRNA\n")
    print "writing combined mRNA-miRNA fold change results in output_combineFC.txt..."

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

"""
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

"""




def main():

    tumorDict, normalDict, total_mirna, allsamplesN, allsamplesT = getFiles_miRNA()
    pval_dict = paired_tTest(total_mirna,allsamplesN,allsamplesT)


    targets_clash = mirnatarget.CLASH()
    targets_mirtar = mirnatarget.MirTarBase()

    outputFCmirna1(targets_clash, pval_dict)

    outputFCmirna2(targets_mirtar,pval_dict)

    mrnadata.main()

    combine_data()

    #correlation()


if __name__ == '__main__':
    main()
