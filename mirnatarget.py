#a script to parse microRNA.org data [and others later on] to find target genes

__author__ = 'sharon'

#processes xls files
import xlrd


#contains miRNA target predictions that are not validated so it's better to use the other DBs
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
            mirna_dict[temp] = (line[2]).lower()
        else:
            mirna_dict[(line[4]).lower()] = (line[2]).lower()
    return mirna_dict


def MirTarBase():
    file = xlrd.open_workbook("//home//sharon//Desktop//miRNA_targets//MirTarBase//hsa_MTI.xls")
    worksheet = file.sheet_by_name('miRTarBase')

    # key: miRNA , keys: [target gene name, target gene ID]
    target_dict = dict()

    for r in range(1,worksheet.nrows):
        #if hit end of excel table
        if worksheet.cell_type(r,1) in (xlrd.XL_CELL_EMPTY,xlrd.XL_CELL_BLANK):
            break
        else:
            mirna = str(worksheet.cell(r,1).value)
            gene_name = str(worksheet.cell(r,3).value)

            #some miRNAs do not have target gene in the excel sheet. only get in the dict the good ones
            gene_ID = worksheet.cell(r,4).value
            if gene_ID != "":
                #print gene_ID
                gene_ID = int(round(float(worksheet.cell(r,4).value)))
                target_dict[(mirna).lower()] = (gene_name.lower(), gene_ID)


    #print target_dict
    return target_dict


def main():

    #microRNA()
    CLASH()
    MirTarBase()



if __name__ == '__main__':
    main()
