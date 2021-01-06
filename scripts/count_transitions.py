#!/usr/bin/python
#-*- coding: utf-8 -*-

#!/usr/bin/python
#-*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import math
import wordanalysis
import os
import re

#Load all files using os and re
##
path = '../Data/Cell_files_data/'
files = os.listdir(path)
filenames = {} #Create a dictionary to save paths for all data cells
for i in files:
    m = re.search(r'.*[^_edited_cells.txt]',i)
    filenames[m.group()] = [i]
my_vocabulary = set()

#Create 
##
for values in filenames.values(): 
    print(values[0], wordanalysis.get_all_words(values[1], 1))

#Quitar los radios y sin que se afecte el orden de las otras células
for values in filenames.values(): 
    values.append([])
    for i in range(0,len(values[1])):
        x = values[1][i].replace('R','')
        if x == '':
            pass
        else:
            values[2].append(x)
