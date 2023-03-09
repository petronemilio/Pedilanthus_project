#!/usr/bin/python
#-*- coding: utf-8 -*-
from collections import defaultdict
import numpy as np
import pandas as pd
import os
import re

path = '../Data/Cell_files_data/ConvergeAssOtherLineage/'
files = os.listdir(path)
filenames = {}
for i in files:
    m = re.search(r'.*[^_edited_cells_NotConverge.txt]',i)
    filenames[m.group()] = [i]
filenames

my_vocabulary = set()
my_lens = {}
for keys, values in filenames.items():
    with open(path+values[0]) as f:
        x = f.readlines()
        filenames.setdefault(keys,[]).append(x)
    values[1] = [x.strip() for x in values[1]]  
    values[1] = [x.upper() for x in values[1]]  
    my_lens[values[0]]=[]
    for s in values[1]:
        my_lens[values[0]].append(len(s))
        for z in s:
            my_vocabulary.add(z)            

mydataframe = pd.DataFrame()
for k in my_lens:
    df = pd.DataFrame()
    df = pd.DataFrame(my_lens[k], columns =['Number of cells'])
    df.insert(0, 'Sample', np.repeat(k, np.shape(df)[0]), True)
    mydataframe = mydataframe.append(df, ignore_index=True)

#Write data fram
mydataframe.to_csv('../Data/cell_lengths_notConverge.csv')
####Do it for cell lineages without R####
my_lens = {}
for values in filenames.values(): 
    values.append([])
    for i in range(0,len(values[1])):
        x = values[1][i].replace('R','')
        if x == '':
            pass
        else:
            values[2].append(x)
    my_lens[values[0]]=[]
    for s in values[2]:
        my_lens[values[0]].append(len(s))
        for z in s:
            my_vocabulary.add(z)            
mydataframe = pd.DataFrame()
for k in my_lens:
    df = pd.DataFrame()
    df = pd.DataFrame(my_lens[k], columns =['Number of cells'])
    df.insert(0, 'Sample', np.repeat(k, np.shape(df)[0]), True)
    mydataframe = mydataframe.append(df, ignore_index=True)

#Write data fram
mydataframe.to_csv('../Data/cell_lengths_withoutR.csv')

