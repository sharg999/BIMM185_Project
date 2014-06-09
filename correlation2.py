

import mirnatarget

import matplotlib.pyplot as plt

__author__ = 'Sharong'

import os
import fnmatch
from scipy import stats


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
                        allsamplesT[line[0]].append(float(0.00000))
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
                        allcountsN[line[0]]= float(0.00000)
                        howmanyN[line[0]] = 1
                        allsamplesN.setdefault(line[0],[])
                        allsamplesN[line[0]].append(float(0.00000))
                    else:
                        allcountsN[line[0]]= float(line[2])
                        howmanyN[line[0]] = 1
                        allsamplesN.setdefault(line[0],[])
                        allsamplesN[line[0]].append(float(line[2]))

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
                if line[0].startswith("gene_id") or line[0].startswith("?"):
                    continue
                else:
                    #id[0] will be gene name, id[1] will be gene id #
                    id = line[0].split('|')
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
                            allcountsT[id[0]]= float(0.00000)
                            howmanyT[id[0]] = 1
                            allsamplesT.setdefault(id[0],[])
                            allsamplesT[id[0]].append(float(0.00000))
                        else:
                            allcountsT[id[0]]= float(norm_count)
                            howmanyT[id[0]] = 1
                            allsamplesT.setdefault(id[0],[])
                            allsamplesT[id[0]].append(float(norm_count))
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
                if line[0].startswith("gene_id") or line[0].startswith("?"):
                    continue
                else:
                    #id[0] will be gene name, id[1] will be gene id #
                    id = line[0].split('|')
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

    return allsamplesT, allsamplesN, gene_count



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
                if l==1:
                    alltargetsT[(tmp).strip()] = v
                elif l>1:
                    for i in range(l):
                        alltargetsT[(tmp[i]).strip()] = v


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
            if l==1:
                alltargetsN[(tmp).strip()] = v
            elif l>1:
                for i in range(l):
                    alltargetsN[(tmp[i]).strip()] = v

    FC_each = dict()
    for (k,v),(k2,v2) in zip(alltargetsT.items(),alltargetsN.items()):
        if k==k2:
            tmp = []
            v = [float(x) for x in v.split(",")]
            v2 = [float(x) for x in v2.split(",")]
            for i in range(0,len(v)):
                tmp.append((float(v[i])+0.1)/((float(v2[i]))+0.1))
            FC_each[k] = tmp


    return FC_each

def pair_mrna(allsamplesT_mrna,allsamplesN_mrna):
    FC_each = dict()
    for (k,v),(k2,v2) in zip(allsamplesT_mrna.items(),allsamplesN_mrna.items()):
        if k==k2:
            tmp = []
            v = [float(x) for x in v]
            v2 = [float(x) for x in v2]
            for i in range(0,len(v)):
                tmp.append((float(v[i])+0.1)/((float(v2[i]))+0.1))
            FC_each[k] = tmp
    return FC_each


def mutual_genes(FC_each_mrna,FC_each_mirna):

    mut_genes = []
    for k in FC_each_mrna.keys():
        for k2 in FC_each_mirna.keys():
            if k==k2:
                #print k
                mut_genes.append(k)
    return mut_genes


def correlation(FC_each_mrna,FC_each_mirna, mut_genes):
    import math
    f = open('//home//sharon//Desktop//TCGA//KICH//correlation.txt', 'w')
    f.write("#The Pearson correlation between mRNA and miRNA log2-fold changes\n")
    f.write("Gene_name\tPearson_coefficient\tP-value\n")
    print "writing correlation results in correlation.txt..."

    d = dict()

    for g in mut_genes:
        x = FC_each_mrna[g]
        y = FC_each_mirna[g]
        tmp, pv = ((stats.spearmanr(y,x)))
        if math.isnan(tmp):
            continue
        else:
            f.write(str(g)+"\t"+str(tmp)+"\t"+str(pv)+"\n")
        if (tmp >0.1) or (tmp <-0.1):
            d[g] = tmp

    l = len(mut_genes)

    plt.title("Pearson Correlation of mRNA vs. miRNA data")
    plt.xlabel("Gene")
    plt.xticks(rotation=90, fontsize = 8, linespacing = float(20.0))

    plt.ylabel("Pearson Correlation Coefficient")

    fig = plt.bar(range(len(d)), d.values(),align='center',linewidth=0.5)
    plt.xticks(range(len(d)),d.keys())
    plt.xlim(xmin=0)
    plt.show()


def main():

    allsamplesN, allsamplesT = getFiles_miRNA()

    targets_clash = mirnatarget.CLASH()
    targets_mirtar = mirnatarget.MirTarBase()
    FC_each_mirna = miRNAtargets(targets_clash,targets_mirtar,allsamplesN,allsamplesT)

    allsamplesT_mrna, allsamplesN_mrna, gene_count = getfiles_mrna()

    FC_each_mrna = pair_mrna(allsamplesT_mrna,allsamplesN_mrna)

    mut_genes = mutual_genes(FC_each_mrna,FC_each_mirna)

    correlation(FC_each_mrna,FC_each_mirna,mut_genes)


if __name__ == '__main__':
    main()
