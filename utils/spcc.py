import pandas as pd
import os
import scipy.stats as stats


"""
The calculation of SPCC
"""
def nonsyn_syn_files():
    filepath = 'D:\Maftools\TMB\WXS_Mutect\Data\Data\Raw Files\\'
    filelist = os.listdir(filepath)
    #这里直接把chr置为index，省去一个for循环
    cytoband = pd.read_csv(r"D:\Maftools\CNV\cytoBand3_exon.txt",sep='\t',index_col=0)
    targetfilepath = r'D:\Maftools\TMB\WXS_Mutect\Data\Data\All\\'
    # targetfilepath = r'D:\Maftools\TMB\WXS_Mutect\Data\Data\Nonsyn_syn Files\\'
    sampledata = pd.read_excel("arm-wide TMB value of each sample.xlsx",sheet_name=None, index_col=0)
    for i in filelist:
        cancername = i.split('.')[0]
        subsamplelist = list(sampledata[cancername].index)
        print(subsamplelist)
        mutation_type = ['Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site','Translation_Start_Site','Silent']
        data = pd.read_csv(filepath+i,sep='\t')[['Hugo_Symbol','Chromosome','Variant_Classification','Start_Position','Tumor_Sample_Barcode']]
        data = data[data['Variant_Classification'].isin(mutation_type)].reset_index(drop=True)
        for j in range(len(data)):
            if(data['Chromosome'][j]!='chrY'):
                if(data['Start_Position'][j]<=cytoband['pend'][data['Chromosome'][j]]):
                    data['Chromosome'][j] = data['Chromosome'][j].replace('chr','')+'p'
                else:
                    data['Chromosome'][j] = data['Chromosome'][j].replace('chr', '') + 'q'
                data['Tumor_Sample_Barcode'][j] = data['Tumor_Sample_Barcode'][j][0:16]
        data = data[data['Tumor_Sample_Barcode'].isin(subsamplelist)].reset_index(drop=True)
        print(cancername)
        print(data)
        data.to_csv(targetfilepath+cancername+'.csv',sep='\t')

def nonsyn_syn_correlation(non_syn_path,dtwpath):
    chr = ["Sample", "1p", "1q", "10p", "10q", "11p", "11q", "12p", "12q",
           "13p", "13q", "14p", "14q", "15p", "15q", "16p","16q", "17p", "17q",
           "18p", "18q", "19p", "19q", "2p", "2q", "20p", "20q", "21p", "21q",
           "22p", "22q", "3p", "3q", "4p", "4q", "5p","5q","6p", "6q", "7p",
           "7q", "8p", "8q", "9p", "9q","Xp","Xq","Total"]
    filelist = os.listdir(non_syn_path)
    res_df = pd.DataFrame()
    count = 0
    for i in filelist:
        nonsyn_syn_data = pd.read_csv(non_syn_path+i,sep='\t')
        dtwdata = pd.read_csv(dtwpath+i,header=None)
        dtwdata.columns = chr
        for j in range(len(dtwdata)):
            dtwdata['Sample'][j] = dtwdata['Sample'][j][0:16]
        commonsample = list(set(nonsyn_syn_data['Sample'])&set(dtwdata['Sample']))
        commonsample.sort()
        nonsyn_syn_data = nonsyn_syn_data[nonsyn_syn_data['Sample'].isin(commonsample)].reset_index(drop=True)
        dtwdata = dtwdata[dtwdata['Sample'].isin(commonsample)].reset_index(drop=True)
        correlation_dict = {}
        correlation_dict['Cancer'] = i.split('.')[0]
        count = 0
        for k in range(1,len(chr)):
            try:
                correlation,p = stats.spearmanr(nonsyn_syn_data[chr[k]],-dtwdata[chr[k]])
                statistic,Wilcoxon_p = stats.wilcoxon(nonsyn_syn_data[chr[k]],-dtwdata[chr[k]],alternative='two-sided')
                if (Wilcoxon_p < 0.05):
                    correlation_dict[k] = round(correlation, 4)
                else:
                    correlation_dict[k] = "NS"
            except ValueError:
                correlation_dict[k] = "NS"
        res_df = res_df.append(pd.DataFrame(correlation_dict,index=[i]))
    res_df.columns = ["Cancer", "1p", "1q", "10p", "10q", "11p", "11q", "12p", "12q",
           "13p", "13q", "14p", "14q", "15p", "15q", "16p","16q", "17p", "17q",
           "18p", "18q", "19p", "19q", "2p", "2q", "20p", "20q", "21p", "21q",
           "22p", "22q", "3p", "3q", "4p", "4q", "5p","5q","6p", "6q", "7p",
           "7q", "8p", "8q", "9p", "9q","Xp","Xq","Total"]
    print(count)
    res_df.to_csv(r'D:\Maftools\TMB\WXS_Mutect\Data\Data\Nonsyn+Syn_wilcoxon_p_DTW_Spearman.csv',sep='\t',index=False)


