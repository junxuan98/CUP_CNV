import pandas as pd
import os

path = r"D:\Maftools\TMB\WXS_Mutect\Data\Data2\\"

"""
the preprocessing of MAF tiles
mark the mutation with p and q
"""
def pqTMB():
    filelist = os.listdir(path)
    cytoband = pd.read_csv(r"D:\Maftools\CNV\cytoBand1.txt", sep='\t', header=None)
    for k in range(0,len(filelist)):
        mutect = pd.read_csv(path+filelist[k],sep='\t')
        mutect = mutect.loc[:, ~mutect.columns.str.contains('^Unnamed')]
        mutect = mutect[['Hugo_Symbol', 'Chromosome', 'Start_Position','End_Position', 'Variant_Classification','Tumor_Sample_Barcode']]
        for i in range(0,len(mutect)):
            for j in range(0,len(cytoband)):
                if (cytoband[0][j] == mutect.loc[i, 'Chromosome']):
                    if (mutect.loc[i, 'Start_Position'] > cytoband[1][j] and mutect.loc[i, 'End_Position'] < cytoband[2][j]):
                        mutect.loc[i, 'Chromosome'] = mutect.loc[i, 'Chromosome'].replace('chr','') + cytoband[3][j][0]
                        break
        mutect.to_csv(path+filelist[k],sep='\t',index = False)

"""
the calculation of exoTMB(nonsyn)
"""
def chrTMB(name):
    filelist = os.listdir('D:\Maftools\CNV\Data\\'+name+'\T_TMB')
    mutect = pd.read_csv(path+name+".maf", sep='\t')
    chr = ["1p","1q","10p","10q","11p","11q","12p","12q","13p","13q","14p","14q","15p","15q","16p","16q","17p","17q",
            "18p","18q","19p","19q","2p","2q","20p","20q","21p","21q","22p","22q","3p","3q","4p","4q","5p","5q",
            "6p","6q","7p","7q","8p","8q","9p","9q"]
    df = pd.DataFrame()
    dict = {}
    for m in range(0, len(filelist)):
        dict['Sample'] = filelist[m][0:16]
        for j in range(0, len(chr)):
            k = 0
            for i in range(0, len(mutect)):
                if (filelist[m][0:16] in mutect['Tumor_Sample_Barcode'][i]):
                    if(mutect['Chromosome'][i].replace('chr','') == chr[j]):
                        if mutect['Variant_Classification'][i] == 'Frame_Shift_Del':
                            k = k + 1
                        if mutect['Variant_Classification'][i] == 'Frame_Shift_Ins':
                            k = k + 1
                        if mutect['Variant_Classification'][i] == 'In_Frame_Del':
                            k = k + 1
                        if mutect['Variant_Classification'][i] == 'In_Frame_Ins':
                            k = k + 1
                        if mutect['Variant_Classification'][i] == 'Missense_Mutation':
                            k = k + 1
                        if mutect['Variant_Classification'][i] == 'Nonsense_Mutation':
                            k = k + 1
                        if mutect['Variant_Classification'][i] == 'Nonstop_Mutation':
                            k = k + 1
                        if mutect['Variant_Classification'][i] == 'Splice_Site':
                            k = k + 1
                        if mutect['Variant_Classification'][i] == 'Translation_Start_Site':
                            k = k + 1
            dict[chr[j]] = '%.2f' % float(k/50)
        df = df.append(pd.DataFrame(dict,index=[m]))
    # df.to_csv('D:\Maftools\TMB\WXS_Mutect\Data\Chr_TMB\\'+name+'_chr_TMB.csv',index=False)
