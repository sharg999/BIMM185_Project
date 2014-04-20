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


def main():

    microRNA()



if __name__ == '__main__':
    main()
