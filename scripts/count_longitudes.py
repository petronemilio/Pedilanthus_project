#!/usr/bin/python
#-*- coding: utf-8 -*-
from collections import defaultdict
import numpy as np
import pandas as pd
import os
import re

path = '../Data/Cell_files_data/'
files = os.listdir(path)
filenames = {}
for i in files:
    m = re.search(r'.*[^_edited_cells.txt]',i)
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
mydataframe.to_csv('../Data/cell_lengths.csv')

