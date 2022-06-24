import pandas as pd
import numpy as np
import glob,os
import json
import re
import time

def newdir(name):
    os.makedirs('D:\\Maftools\\CNV\\Data\\' + name + "\\N")
    os.makedirs('D:\\Maftools\\CNV\\Data\\' + name + "\\T")
    os.makedirs('D:\\Maftools\\CNV\\Data\\' + name + "\\N_TMB")
    os.makedirs('D:\\Maftools\\CNV\\Data\\' + name + "\\T_TMB")
    os.makedirs('D:\\Maftools\\CNV\\Data\\' + name + "\\Data")
    os.makedirs('D:\\Maftools\\CNV\\Data\\' + name + "\\FinalData")

def readjson(filepath):
    file = open(filepath,'r',encoding = 'utf-8' )
    data = json.load(file)
    c = []
    for i in range(0,len(data)):
        c.append(data[i]['associated_entities'][0]['entity_submitter_id']+" "+data[i]['file_name'])
    return c

def renamefiles(name):
    filelist = os.listdir('D:\\Maftools\\CNV\\Data\\' + name + '\\Raw_data\\')
    metadata = readjson("D:\\Maftools\\CNV\\Data\\" + name + "\\metadata.json")
    sheet = pd.read_csv(r'D:\\Maftools\\CNV\\Data\\' + name + '\\sheet.tsv', sep='\t')
    k = 0
    for i in range(0, len(metadata)):
        if (metadata[i].split(' ')[1] == sheet.loc[i, 'File Name']):
            k = k + 1
    if (k == len(metadata)):
        print('metadata 和 sheet 的文件名一一对应')
        for i in range(0, len(filelist)):
            for j in range(0, len(metadata)):
                if filelist[i] == metadata[j].split(' ')[1]:
                    if (sheet.loc[j, 'Sample Type'] == 'Primary Tumor' or sheet.loc[j, 'Sample Type'] == 'Solid Tissue Normal'
                    or sheet.loc[j, 'Sample Type'] == 'Blood Derived Normal'):
                        olddir = 'D:\\Maftools\\CNV\\Data\\' + name + '\\Raw_data\\' + filelist[i]
                        newdir = 'D:\\Maftools\\CNV\\Data\\' + name + '\\Raw_data\\' + metadata[j].split(' ')[0] + '-' + \
                                 sheet.loc[j, 'Sample Type'] + '.txt'
                        os.rename(olddir, newdir)
                    else:
                        os.remove('D:\\Maftools\\CNV\\Data\\' + name + '\\Raw_data\\' + filelist[i])
        print('rename success')
    else:
        print('metadata can not meet sheet ')


