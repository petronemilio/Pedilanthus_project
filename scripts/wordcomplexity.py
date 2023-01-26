#!/usr/bin/python
#-*- coding: utf-8 -*-

#import matplotlib.pyplot as plt
#import numpy as np
import math
import wordanalysis
import os
import re
import pandas as pd
import numpy as np
import itertools
import csv

#Load all files using os and re
##
path = '../Data/Cell_files_data/ConvergeAssOtherLineage/'
files = os.listdir(path)
filenames = {} #Create a dictionary to save paths for all data cells
for i in files:
    m = re.search(r'.*[^_edited_cells_NotConverge.txt]',i)
    filenames[m.group()] = [i]

#Loop to load the files of cells and append them in filenames
my_vocabulary = set()
for keys, values in filenames.items():
    with open(path+values[0]) as f:
        x = f.readlines()
        filenames.setdefault(keys,[]).append(x)
    values[1] = [x.strip() for x in values[1]]  #remove the /n 
    values[1] = [x.upper() for x in values[1]]  #put all letters in upper
    for s in values[1]:
        for z in s:
            my_vocabulary.add(z)

#Remove R
for values in filenames.values(): 
    values.append([])
    for i in range(0,len(values[1])):
        x = values[1][i].replace('R','')
        if x == '':
            pass
        else:
            values[2].append(x)


lempelziv_log = {} #Empty dictionary to add values into
for values in filenames.values():
    lempelziv_log[values[0]]=[]
    for z in range(len(values[2])):
        x = wordanalysis.lzcomprdist(values[2][z])
        lempelziv_log[values[0]].append(x) 

mydataframe = pd.DataFrame()
for k,v in lempelziv_log.items():
    lst2 = list(np.repeat(k,len(v)))
    df = pd.DataFrame(list(zip(v, lst2)), columns =['Value', 'Name'])
    #df.insert(0, 'Sample', np.repeat(k, np.shape(df)[0]), True)
    mydataframe = mydataframe.append(df, ignore_index=True)

mydataframe.to_csv('../Data/lemplzivbyfile.csv')


######Calc. shannon entropy of individual cell files

shannon_entropy_files = {} #Empty dictionary to add values into
for values in filenames.values():
    shannon_entropy_files[values[0]]=[]
    for z in range(len(values[2])):
        x = wordanalysis.estimate_shannon_entropy(values[2][z])
        shannon_entropy_files[values[0]].append(x)

mydataframe = pd.DataFrame()
for k,v in shannon_entropy_files.items():
    lst2 = list(np.repeat(k,len(v)))
    df = pd.DataFrame(list(zip(v, lst2)), columns =['Value', 'Name'])
    #df.insert(0, 'Sample', np.repeat(k, np.shape(df)[0]), True)
    mydataframe = mydataframe.append(df, ignore_index=True)

mydataframe.to_csv('../Data/shannonentropy.csv')

#######Calc. shannon entropy of cell files considering window 
shannon_by_window={}
for values in filenames.values():
    shannon_by_window[values[0]]={}
    x=1
    for z in values[2]:
        shannon_by_window[values[0]][x] = wordanalysis.estimate_shannon_entropy_window(z)
        x+=1

#Sigue exportar a data frmae o algo parecido con pandas
mydataframe = pd.DataFrame()
for k,v in shannon_by_window.items():
    for z,y in v.items():
        lst2 = list(np.repeat(z,len(y)))
        lst3 = list(np.repeat(k,len(y)))
        df = pd.DataFrame(list(zip(lst3,y, lst2)), columns =['Sample','Value', 'Lineage'])
    #df.insert(0, 'Sample', np.repeat(k, np.shape(df)[0]), True)
        mydataframe = mydataframe.append(df, ignore_index=True)


mydataframe.to_csv('../Data/shannonentropy_window.csv')


