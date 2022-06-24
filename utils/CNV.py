import pandas as pd
import numpy as np
import os
import re
import time

"""
The preprocessing of CNA data
Arrange the positive and negative samples raw data to data matrix respectively
The data matrix row->segment col->samples
"""
def cnv(cancername,type):
    gene = pd.read_csv("D:\\Maftools\\CNV\\b.txt",sep='\t')
    time_start=time.time()
    filepath = 'D:\\Maftools\\CNV\\Data\\' + cancername + '\\' + type + "_TMB\\"
    filelist = os.listdir(filepath)
    list = []
    initialization = gene.values.tolist()
    for i in range(0,1):
        filedata = pd.read_csv(filepath+filelist[i],sep='\t')
        fl = filedata.values.tolist()
        for j in initialization:
            for k in fl:
                if(str(j[0])== re.findall("\d+",k[0])[0]):
                    if(j[2]<= k[1]):
                        j.append(np.nan)
                        break
                    if(j[1]>= k[1] and j[2]<= k[2]):
                        j.append(k[3])
                        break
                    if(j[1]>=k[2]):
                        continue
                    else:
                        j.append(np.nan)
                        break
            if(str(j[0])== re.findall("\d+",k[0])[0]):
                j[0] = str(j[0])+re.findall(r'[a-z]',k[0])[0]
            list.append(j)
    for i in range(1,len(filelist)):
        filedata = pd.read_csv(r"D:\\Maftools\\CNV\\Data\\"+cancername+"\\"+type+"_TMB\\" + filelist[i], sep='\t')
        fl = filedata.values.tolist()
        for j in range(0, len(initialization)):
            for k in range(0, len(fl)):
                if(initialization[j][0] == fl[k][0]):
                    if (initialization[j][2] <= fl[k][1]):
                        list[j].append(np.nan)
                        break
                    if (initialization[j][1] >= fl[k][1] and initialization[j][2] <= fl[k][2]):
                        list[j].append(fl[k][3])
                        break
                    if (initialization[j][1]>=fl[k][2]):
                        continue
                    else:
                        list[j].append(np.nan)
                        break
    data = pd.DataFrame(data=list)
    data = data.dropna(thresh=data.shape[1]-5)
    data.to_csv('D:\\Maftools\\CNV\\Data\\'+cancername+'\\FinalData\\'+cancername+' '+type+'.csv',sep='\t',index = False)
    time_end=time.time()
    print('totally cost',time_end-time_start)








