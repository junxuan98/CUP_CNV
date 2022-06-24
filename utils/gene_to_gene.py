import pandas as pd
import os


#the union set of tsample and its corresponding nsample
def gene_to_gene(name):
    nfilepath = "D:\\Maftools\\CNV\\Data\\" + name + "\\N_TMB\\"
    tfilepath = "D:\\Maftools\\CNV\\Data\\" + name + "\\T_TMB\\"
    nfilelist = os.listdir(nfilepath)
    tfilelist = os.listdir(tfilepath)
    for i in range(0,len(nfilelist)):
        ndf = pd.DataFrame()
        tdf = pd.DataFrame()
        ndict = {}
        tdict = {}
        ndata = pd.read_csv(nfilepath+nfilelist[i],sep='\t')
        ndata = (ndata[ndata.Chromosome == '3q']).reset_index(drop=True)
        tdata = pd.read_csv(tfilepath + tfilelist[i], sep='\t')
        tdata = (tdata[tdata.Chromosome == '3q']).reset_index(drop=True)
        for j in range(0,len(ndata)):
            for k in range(0,len(tdata)):
                # t   |_______|
                # n |____________|
                if(tdata.Start[k]>=ndata.Start[j] and tdata.End[k]<=ndata.End[j]):
                    ndict['Chromosome'] = ndata.Chromosome[j]
                    tdict['Chromosome'] = tdata.Chromosome[k]
                    ndict['Start'] = tdata.Start[k]
                    tdict['Start'] = tdata.Start[k]
                    ndict['End'] = tdata.End[k]
                    tdict['End'] = tdata.End[k]
                    ndict['Segment_Mean'] = ndata.Segment_Mean[j]
                    tdict['Segment_Mean'] = tdata.Segment_Mean[k]
                    if ndict['Start'] != ndict['End']:
                        ndf = ndf.append(pd.DataFrame(ndict,index=[k]))
                        tdf = tdf.append(pd.DataFrame(tdict, index=[k]))
                    continue
                # t           |_______|
                # n |____________|
                elif ((tdata.Start[k]>=ndata.Start[j] and tdata.End[k]>=ndata.End[j]) and ndata.End[j]>=tdata.Start[k]):
                    ndict['Chromosome'] = ndata.Chromosome[j]
                    tdict['Chromosome'] = tdata.Chromosome[k]
                    ndict['Start'] = tdata.Start[k]
                    tdict['Start'] = tdata.Start[k]
                    ndict['End'] = ndata.End[j]
                    tdict['End'] = ndata.End[j]
                    ndict['Segment_Mean'] = ndata.Segment_Mean[j]
                    tdict['Segment_Mean'] = tdata.Segment_Mean[k]
                    if ndict['Start'] != ndict['End']:
                        ndf = ndf.append(pd.DataFrame(ndict, index=[k]))
                        tdf = tdf.append(pd.DataFrame(tdict, index=[k]))
                    break
                # t           |_______|
                # n |____________|  |____________|
                elif ((tdata.Start[k]<=ndata.Start[j] and tdata.End[k] >=ndata.Start[j]) and tdata.End[k]<=ndata.End[j]):
                    ndict['Chromosome'] = ndata.Chromosome[j]
                    tdict['Chromosome'] = tdata.Chromosome[k]
                    ndict['Start'] = ndata.Start[j]
                    tdict['Start'] = ndata.Start[j]
                    ndict['End'] = tdata.End[k]
                    tdict['End'] = tdata.End[k]
                    ndict['Segment_Mean'] = ndata.Segment_Mean[j]
                    tdict['Segment_Mean'] = tdata.Segment_Mean[k]
                    if ndict['Start'] != ndict['End']:
                        ndf = ndf.append(pd.DataFrame(ndict, index=[k]))
                        tdf = tdf.append(pd.DataFrame(tdict, index=[k]))
                    continue
                # t           |_____________|
                # n |____________|  |_____|
                elif (tdata.Start[k] <= ndata.Start[j] and tdata.End[k] >= ndata.End[j]) :
                    ndict['Chromosome'] = ndata.Chromosome[j]
                    tdict['Chromosome'] = tdata.Chromosome[k]
                    ndict['Start'] = ndata.Start[j]
                    tdict['Start'] = ndata.Start[j]
                    ndict['End'] = ndata.End[j]
                    tdict['End'] = ndata.End[j]
                    ndict['Segment_Mean'] = ndata.Segment_Mean[j]
                    tdict['Segment_Mean'] = tdata.Segment_Mean[k]
                    if ndict['Start'] != ndict['End']:
                        ndf = ndf.append(pd.DataFrame(ndict, index=[k]))
                        tdf = tdf.append(pd.DataFrame(tdict, index=[k]))
                    break
                elif (tdata.Start[k]>=ndata.End[j]):
                    break
                elif (tdata.End[k]<=ndata.Start[j]):
                    continue

        ndf.to_csv('D:\Maftools\CNV\Data\LUSC\\3_N\\'+nfilelist[i],index=False)
        tdf.to_csv('D:\Maftools\CNV\Data\LUSC\\3_T\\'+tfilelist[i],index=False)

if __name__ == '__main__':
    gene_to_gene("LUSC")



