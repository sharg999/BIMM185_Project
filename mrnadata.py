__author__ = 'sharon'

import os
import fnmatch
from scipy import stats
import numpy as np

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
        #parse files for tumor samples
        if fnmatch.fnmatch(file, '*.rsem.genes.normalized_results') and (tpat==filemap_dict[file][13:15]):

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
                if line[0].startswith("gene_id"):
                    continue
                else:
                    #id[0] will be gene name, id[1] will be gene id #
                    id = line[0].split('|')
                    id[1]= id[1][1:]
                    norm_count = line[1]

                    #otherwise, just add to dictionary and keep track of sample counts
                    if id[1] in allcountsT:
                        #if norm_count< .125 set to 0
                        if float(norm_count) <0.125:
                            allcountsT[(id[0],id[1])]+= float(0.00000)
                            howmanyT[id[1]]+=1
                            allsamplesT[(id[0],id[1])].append(float(0.00000))
                        else:
                            allcountsT[(id[0],id[1])]+= float(norm_count)
                            howmanyT[(id[0],id[1])]+=1
                            allsamplesT[(id[0],id[1])].append(float(norm_count))
                    else:
                        if float(norm_count) < 0.125:
                            allcountsT[(id[0],id[1])]= float(0.0000)
                            howmanyT[(id[0],id[1])] = 1
                            allsamplesT.setdefault((id[0],id[1]),[])
                            allsamplesT[(id[0],id[1])].append(float(0.0000))
                        else:
                            allcountsT[(id[0],id[1])]= float(norm_count)
                            howmanyT[(id[0],id[1])] = 1
                            allsamplesT.setdefault((id[0],id[1]),[])
                            allsamplesT[(id[0],id[1])].append(float(norm_count))
        elif fnmatch.fnmatch(file, '*.rsem.genes.normalized_results') and ('11'==filemap_dict[file][13:15]):

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
                if line[0].startswith("gene_id"):
                    continue
                else:
                    #id[0] will be gene name, id[1] will be gene id #
                    id = line[0].split('|')
                    id[1]= id[1][1:]
                    norm_count = line[1]

                    #otherwise, just add to dictionary and keep track of sample counts
                    if id[1] in allcountsN:
                        #if norm_count< .125 set to 0
                        if float(norm_count) <0.125:
                            allcountsN[(id[0],id[1])]+= float(0.00000)
                            howmanyN[id[1]]+=1
                            allsamplesN[(id[0],id[1])].append(float(0.00000))
                        else:
                            allcountsN[(id[0],id[1])]+= float(norm_count)
                            howmanyN[(id[0],id[1])]+=1
                            allsamplesN[(id[0],id[1])].append(float(norm_count))
                    else:
                        if float(norm_count) < 0.125:
                            allcountsN[(id[0],id[1])]= float(0.0000)
                            howmanyN[(id[0],id[1])] = 1
                            allsamplesN.setdefault((id[0],id[1]),[])
                            allsamplesN[(id[0],id[1])].append(float(0.0000))
                        else:
                            allcountsN[(id[0],id[1])]= float(norm_count)
                            howmanyN[(id[0],id[1])] = 1
                            allsamplesN.setdefault((id[0],id[1]),[])
                            allsamplesN[(id[0],id[1])].append(float(norm_count))

    gene_count = len(allsamplesT)

    return allsamplesT, allsamplesN, gene_count


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
            foldchange_dict[(k[0],k[1])] = float(0)
        else:
            foldchange_dict[(k[0],k[1])] = np.log2((float(t)+0.1) / ((float(n)+0.1)))
        #if foldchange_dict[(k[0],k[1])] >=1.0 or foldchange_dict[(k[0],k[1])]<=-1.0:
        result[k] = foldchange_dict[k]

    return result


def paired_tTest(total_genes, allsamplesN, allsamplesT):

    #output of all differentially expressed genes
    #output = open('//home//sharon//Desktop//TCGA//BRCA//output_mRNA_genes.txt', 'w')
    #output = open('//home//sharon//Desktop//TCGA//LIHC//output_mRNA_genes.txt', 'w')
    #output = open('//home//sharon//Desktop//TCGA//LUAD//output_mRNA_genes.txt', 'w')
    #output = open('//home//sharon//Desktop//TCGA//ESCA//output_mRNA_genes.txt', 'w')
    #output = open('//home//sharon//Desktop//TCGA//HNSC//output_mRNA_genes.txt', 'w')
    output = open('//home//sharon//Desktop//TCGA//KICH//output_mRNA_genes.txt', 'w')


    output.writelines("#differentially expressed genes from mRNA data. FC threshold -+1, p_val<0.05 \n")
    output.writelines("gene_name"+"\t"+"gene_id"+"\t"+"log2_FC"+"\t"+"p_val"+"\n")

    #stores all the significant p-vals <0.05
    pval_dict = dict()

    FC = foldChange(allsamplesN,allsamplesT)

    n = total_genes
    significant_count = 0
    for k,v in allsamplesT.items():
        if k in allsamplesN:
            normal = allsamplesN[k]
            tumor = allsamplesT[k]

            paired_sample = stats.ttest_rel(normal,tumor)

            #signficant if p<0.05 and FC = +-1
            #if paired_sample[1] < 0.05 and k in FC:
            if (FC[k]>=1.0 or FC[k]<=-1.0 ) and paired_sample[1] <0.05:
                    #print("paired sample:", paired_sample)
                    significant_count +=1
                    output.write(k[0].lower()+"\t"+k[1] + "\t"+ str(FC[k])+"\t"+str(paired_sample[1])+"\n")
    print "Number of significant differentially expressed miRNA targets: ", significant_count
    output.close()
    return pval_dict


def main():
    allsamplesT, allsamplesN, gene_count = getfiles_mrna()
    pval = paired_tTest(gene_count,allsamplesN,allsamplesT)


if __name__ == '__main__':
    main()
