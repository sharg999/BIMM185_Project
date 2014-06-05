__author__ = 'sharon'


#parse the paired samples from both miRNA and mRNA for pearson correlation


import os
import fnmatch
from scipy import stats
import numpy as np
import mirnatarget

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
            # {(barcode) : file name}
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
        #print file
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

                    #normalized count of gene
                    norm_count = line[1]

                    #patient
                    patient = filemap_dict[file][8:12]


                    #otherwise, just add to dictionary and keep track of sample counts
                    if id[1] in allcountsT:
                        #if norm_count< .125 set to 0
                        if float(norm_count) <0.125:
                            allcountsT[(id[0],id[1])]+= float(0.00000)
                            howmanyT[id[1]]+=1
                            allsamplesT[(patient, id[0],id[1])].append(float(0.00000))
                        else:
                            allcountsT[(id[0],id[1])]+= float(norm_count)
                            howmanyT[(id[0],id[1])]+=1
                            allsamplesT[(patient, id[0],id[1])].append(float(norm_count))
                    else:
                        if float(norm_count) < 0.125:
                            allcountsT[(id[0],id[1])]= float(norm_count)
                            howmanyT[(id[0],id[1])] = 1
                            allsamplesT.setdefault((patient,id[0],id[1]),[])
                            allsamplesT[(patient,id[0],id[1])].append(float(norm_count))
                        else:
                            allcountsT[(id[0],id[1])]= float(norm_count)
                            howmanyT[(id[0],id[1])] = 1
                            allsamplesT.setdefault((patient,id[0],id[1]),[])
                            allsamplesT[(patient,id[0],id[1])].append(float(norm_count))
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
                if line[0].startswith("gene_id"):
                    continue

                else:
                    #id[0] will be gene name, id[1] will be gene id #
                    id = line[0].split('|')
                    id[1]= id[1][1:]


                    #normalized gene count
                    norm_count = line[1]

                    #patient
                    patient = filemap_dict[file][8:12]


                    #otherwise, just add to dictionary and keep track of sample counts
                    if id[1] in allcountsN:
                        #if norm_count< .125 set to 0
                        if float(norm_count) <0.125:
                            allcountsN[(id[0],id[1])]+= float(0.00000)
                            howmanyN[id[1]]+=1
                            allsamplesN[(patient,id[0],id[1])].append(float(0.00000))
                        else:
                            allcountsN[(id[0],id[1])]+= float(norm_count)
                            howmanyN[(id[0],id[1])]+=1
                            allsamplesN[(patient,id[0],id[1])].append(float(norm_count))
                    else:
                        if float(norm_count) < 0.125:
                            allcountsN[(id[0],id[1])]= float(norm_count)
                            howmanyN[(id[0],id[1])] = 1
                            allsamplesN.setdefault((patient,id[0],id[1]),[])
                            allsamplesN[(patient,id[0],id[1])].append(float(norm_count))
                        else:
                            allcountsN[(id[0],id[1])]= float(norm_count)
                            howmanyN[(id[0],id[1])] = 1
                            allsamplesN.setdefault((patient,id[0],id[1]),[])
                            allsamplesN[(patient,id[0],id[1])].append(float(norm_count))


    gene_count = len(allsamplesT)


    return allsamplesT, allsamplesN


##########

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

    #allcountsN1 = dict()
    #allcountsT1 = dict()

    #counts how many miRNA data I have per each
    howmanyN = dict()
    howmanyT = dict()

    #01 at a specific location means tumor sample
    tpat = "01"

    #dictionary for normal and tumor samples-- keeps all values
    allsamplesN = dict()
    allsamplesT = dict()

    for file in os.listdir(path):

        patient = file[8:12]
        #parse files for tumor samples
        if fnmatch.fnmatch(file, '*.mirna.quantification.txt') and (tpat==file[13:15]):
            #print(file[13:15])

            #f = open('//home//sharon//Desktop//TCGA//BRCA//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//LIHC//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//LUAD//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//ESCA//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//HNSC//miRNA//%s' % file, 'r')
            f = open('//home//sharon//Desktop//TCGA//KICH//miRNA//%s' % file, 'r')

            #patient = file[8:12]




            lines = f.readlines()
            #print(lines)
            for l in lines:
                #print allcountsT


                line = l.split('\t')
                line= list(line)
                #if it's the file title just ignore
                if line[0].startswith("miRNA_ID"):
                    continue
                #otherwise, just add to dictionary and keep track of sample counts

                elif (patient,line[0]) in allcountsT:

                    #if RPM< .125 set to 0
                    if float(line[2]) <0.125:
                        allcountsT[(patient,line[0])]+= float(0.00000)
                        ###allcountsT1[patient]+= float(0.00000)

                        howmanyT[line[0]]+=1
                        allsamplesT[(patient,line[0])].append(float(0.00000))
                    else:
                        allcountsT[(patient,line[0])]+= float(line[2])
                        ###allcountsT1[patient]+= float(line[2])

                        howmanyT[line[0]]+=1

                        allsamplesT[(patient,line[0])].append(float(line[2]))
                else:
                    if float(line[2]) < 0.125:
                        allcountsT[(patient,line[0])]= float(0.0000)
                        ###allcountsT1[patient]= float(0.0000)
                        howmanyT[line[0]] = 1
                        allsamplesT.setdefault((patient,line[0]),[])
                        allsamplesT[(patient,line[0])].append(float(line[2]))
                    else:
                        allcountsT[(patient,line[0])]= float(line[2])
                        ###allcountsT1[patient]= float(line[2])
                        howmanyT[line[0]] = 1
                        allsamplesT.setdefault((patient,line[0]),[])
                        allsamplesT[(patient,line[0])].append(float(line[2]))

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


                #print allcountsN1
                #if it's the file title just ignore
                if line[0].startswith("miRNA_ID"):
                    continue
                #otherwise, just add to dictionary and keep track of sample counts
                elif (patient,line[0]) in allcountsN:
                    #if RPM < 0.125 set to 0
                    if float(line[2]) < 0.125:
                        allcountsN[(patient,line[0])]+= float(0.00000)
                        ###allcountsN1[patient]+= float(0.00000)
                        howmanyN[line[0]]+=1
                        allsamplesN[line[0]].append(float(0.00000))
                    else:
                       allcountsN[(patient,line[0])]+= float(line[2])
                       ###allcountsN1[patient]+= float(line[2])
                       howmanyN[line[0]]+=1
                       allsamplesN[line[0]].append(float(line[2]))
                else:
                    if float(line[2]) < 0.125:
                        allcountsN[(patient,line[0])]= float(line[2])
                        ###allcountsN1[patient]= float(line[2])
                        howmanyN[line[0]] = 1
                        allsamplesN.setdefault(line[0],[])
                        allsamplesN[line[0]].append(float(line[2]))
                    else:
                        allcountsN[(patient,line[0])]= float(line[2])
                        ###allcountsN1[patient]= float(line[2])
                        howmanyN[line[0]] = 1
                        allsamplesN.setdefault(line[0],[])
                        allsamplesN[line[0]].append(float(line[2]))



    #allcounts have the total RPM for that miRNA divided by the number of total samples
    #allsamples includes all the RPM values collected per each miRNA
    return  allsamplesN,allsamplesT



"""
def getfiles_mirna():

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
    allsampleN = dict()
    allsampleT = dict()

    for file in os.listdir(path):
        #parse files for tumor samples
        if fnmatch.fnmatch(file, '*.mirna.quantification.txt') and (tpat==file[13:15]):

            #f = open('//home//sharon//Desktop//TCGA//BRCA//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//LIHC//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//LUAD//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//ESCA//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//HNSC//miRNA//%s' % file, 'r')
            f = open('//home//sharon//Desktop//TCGA//KICH//miRNA//%s' % file, 'r')

            #patient
            patient = file[8:12]
            print(patient)


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
                        allsampleT[(patient,line[0])].append(float(0.00000))
                    else:
                        allcountsT[line[0]]+= float(line[2])
                        howmanyT[line[0]]+=1
                        allsampleT[(patient,line[0])].append(float(line[2]))
                else:
                    if float(line[2]) < 0.125:
                        allcountsT[line[0]]= float(line[2])
                        howmanyT[line[0]] = 1
                        allsampleT.setdefault((patient,line[0]),[])
                        allsampleT[(patient,line[0])].append(float(line[2]))
                    else:
                        allcountsT[line[0]]= float(line[2])
                        howmanyT[line[0]] = 1
                        allsampleT.setdefault((patient,line[0]),[])
                        allsampleT[(patient,line[0])].append(float(line[2]))

        #parse the files for normal samples
        elif fnmatch.fnmatch(file, '*.mirna.quantification.txt') and (tpat != file[13:15]):


            #f = open('//home//sharon//Desktop//TCGA//BRCA//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//LIHC//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//LUAD//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//ESCA//miRNA//%s' % file, 'r')
            #f = open('//home//sharon//Desktop//TCGA//HNSC//miRNA//%s' % file, 'r')
            f = open('//home//sharon//Desktop//TCGA//KICH//miRNA//%s' % file, 'r')

            #patient
            patient = file[8:12]
            print(patient)

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
                        allsampleN[(patient,line[0])].append(float(0.00000))
                    else:
                       allcountsN[line[0]]+= float(line[2])
                       howmanyN[line[0]]+=1
                       allsampleN[(patient,line[0])].append(float(line[2]))
                else:
                    if float(line[2]) < 0.125:
                        allcountsN[line[0]]= float(line[2])
                        howmanyN[line[0]] = 1
                        allsampleN.setdefault((patient,line[0]),[])
                        allsampleN[(patient,line[0])].append(float(line[2]))
                    else:
                        allcountsN[line[0]]= float(line[2])
                        howmanyN[line[0]] = 1
                        allsampleN.setdefault((patient,line[0]),[])
                        allsampleN[(patient,line[0])].append(float(line[2]))


    mirna_count = len(allcountsT)


    #allcounts have the total RPM for that miRNA divided by the number of total samples
    #allsamples includes all the RPM values collected per each miRNA
    return allsampleN,allsampleT
"""


def getFC_mRNA(allsamplesT, allsamplesN):

    #get the FC for each gene
    geneFC = dict()

    #print allsamplesN
    #print allsamplesT

    for (k,v), (k2,v2) in zip(allsamplesT.items(),allsamplesN.items()):
        if k[0] == k2[0]:

            gene_name = k[1]
            #print str(v2)[1:-1]
            count_normal = float(str(v2)[1:-1])
            count_tumor = float(str(v)[1:-1])
            geneFC[gene_name] = np.log2(float((count_tumor+0.1) / (count_normal+0.1)))
    print geneFC

    """
    #pearson:
    result_pearson =  stats.pearsonr(x,y)
    print(result_pearson)
    """

"""
    for k,v in allsamplesT.items():
        patient = k[0]
        gene_num = k[2]
        for k
        if patient in allsamplesN.keys():
            val = allsamplesN.get(patient)

            count_normal = float(allsamplesN[gene_num])
            count_tumor = float(allsamplesT[gene_num])
            #calculate the FC for that gene
            #WHAT ABOUT DIVING BY ZERO?!?!?!
            geneFC[k[1]] = np.log2(float((count_tumor+0.1) / (count_normal+0.1)))


    print geneFC
"""


##def miRNAtargets(targets_clash,targets_mirtar,allsamplesN_mirna,allsamplesT_mirna):
def miRNAtargets(targets_clash,targets_mirtar):

    #process all miRNA targets to once dictionary
    genetargets = dict([],)
    for (k,v),(k2,v2) in zip(targets_clash.items(),targets_mirtar.items()):
        print genetargets
        if k not in genetargets:
            genetargets.setdefault(k,[])
            genetargets[k].append(v)

        else:
            if v in genetargets[k]:
                continue
            else:
                genetargets[k].append(v)

        if k2 not in genetargets:
            genetargets.setdefault(k2,[])
            genetargets[k2] = str(v2[0])
        else:
            print k2,v2[0]

            if v2[0] in genetargets[k2]:
                continue
            else:
                genetargets[k2].append(str(v2[0]))



    print genetargets
    #tumor targets



def main():
    ###allsamplesT_mrna, allsamplesN_mrna = getfiles_mrna()

    ###allsamplesN_mirna, allsamplesT_mirna = getFiles_miRNA()
    #allsamplesN_mirna,allsamplesT_mirna = getfiles_mirna()
    targets_clash = mirnatarget.CLASH()
    targets_mirtar = mirnatarget.MirTarBase()

    #print targets_mirtar
    #print targets_clash

    ##miRNAtargets(targets_clash,targets_mirtar,allsamplesN_mirna,allsamplesT_mirna)
    miRNAtargets(targets_clash,targets_mirtar)

    #print allsamplesT_mirna



    getFC_mRNA(allsamplesT_mrna,allsamplesN_mrna)



if __name__ == '__main__':
    main()
