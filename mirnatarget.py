#a script to parse microRNA.org data [and others later on] to find target genes

__author__ = 'sharon'


def microRNA():
    file = open("//home//sharon//Desktop//miRNA_targets//microrna.org//hg19_predictions_S_C_aug2010.txt", 'r')
    title = file.readline().strip()
    mirna_dict = dict()
    for line in file:
        line = line.strip()
        line = line.split("\t")
        mirna_dict[((line[1]).lower(),(line[2]).lower(),(line[3]).lower())] = float(line[18])
    return mirna_dict


def CLASH():
    file = open("//home//sharon//Desktop//miRNA_targets//CLASH//target_miRSNP_human_CLASH.txt", 'r')
    title = file.readline().strip()
    mirna_dict = dict()
    for line in file:
        line = line.strip()
        line = line.split("\t")
        ## delete the astreisk from the miRNA name
        if line[4][-1] == '*':
            temp = (line[4][:-1]).lower()
            #print((line[4]).lower())
            #print (line[4][:-1]).lower()
            mirna_dict[temp] = (line[2]).lower()
        else:
            mirna_dict[(line[4]).lower()] = (line[2]).lower()
    return mirna_dict


def main():

    microRNA()
    CLASH()



if __name__ == '__main__':
    main()
