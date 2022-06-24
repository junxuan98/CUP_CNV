import pandas as pd
import os

"""
preprocessing : mark the raw data with p and q
"""

def pq(name):
    cytoband = pd.read_csv(r"D:\Maftools\CNV\cytoBand1.txt",sep='\t',header=None)
    filelist = os.listdir("D:\\Maftools\\CNV\\Data\\"+name+"\\Raw_data")
    for k in range(0,len(filelist)):
        data = pd.read_csv('D:\\Maftools\\CNV\\Data\\'+name+'\\Raw_data\\'+filelist[k],sep='\t')
        data = data.drop(columns=['GDC_Aliquot','Num_Probes'])
        data = data[~data['Chromosome'].str.contains("Y")]
        for i in range(0,len(data)):
            for j in range(0,len(cytoband)):
                if (cytoband[0][j].replace('chr','') == data.loc[i,'Chromosome']):
                    if(data.loc[i,'Start']>cytoband[1][j] and data.loc[i,'End']<cytoband[2][j]):
                        data.loc[i,'Chromosome'] = data.loc[i,'Chromosome']+cytoband[3][j][0]
                        break
                    if((data.loc[i,'Start']<cytoband[2][j]) and (data.loc[i,'End']>cytoband[1][j+1] and data.loc[i,'End']<cytoband[2][j+1])):
                        row = {'Chromosome':[data.loc[i,'Chromosome']+cytoband[3][j][0],data.loc[i,'Chromosome']+cytoband[3][j+1][0]],
                               'Start':[data.loc[i,'Start'],cytoband[2][j]],
                               'End':[cytoband[1][j+1],data.loc[i,'End']],
                               'Segment_Mean':[data.loc[i,'Segment_Mean'],data.loc[i,'Segment_Mean']]}
                        df = pd.DataFrame(data = row)
                        data = data.append(df,ignore_index=True)
                        break
        for i in range(0,len(data)):
            if(data.loc[i,'Chromosome'].isnumeric()):
                data =data.drop(i)
        data = data.sort_values(by=['Chromosome','Start','End'])


