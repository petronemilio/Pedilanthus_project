#!/usr/bin/python
#-*- coding: utf-8 -*-

#import matplotlib.pyplot as plt
#import numpy as np
import math
import wordanalysis
import os
import re
import pandas as pd
import itertools
import csv
#Load all files using os and re
##
path = '../Data/Cell_files_data/'
files = os.listdir(path)
filenames = {} #Create a dictionary to save paths for all data cells
for i in files:
    m = re.search(r'.*[^_edited_cells.txt]',i)
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
            
#Now dictoniary is ready for analysis.
#Quitar los radios y sin que se afecte el orden de las otras células
for values in filenames.values(): 
    values.append([])
    for i in range(0,len(values[1])):
        x = values[1][i].replace('R','')
        if x == '':
            pass
        else:
            values[2].append(x)

#################  Create a dictionary to save the word counting #########
word_diversity = {} #Empty dictionary to add values into
for i in range(1,11):
    word_diversity[i]={}
    for values in filenames.values(): 
        word_diversity[i][values[0]] = wordanalysis.get_all_words(values[2], i)

###Write matrix output from df
val = 1
for x in word_diversity.values():
    wordcount = pd.DataFrame.from_dict(x, orient='index')
    wordcount.to_csv('../Data/wordcounts'+ str(val)+'.csv', index = True)    
    val += 1

################################# Calc Euclidean distance   #################
euc_list = {}
for k, v in word_diversity.items():
    header = ["File1","File2", "Euc_dist"]
    euc_list[k]={}
    my_list=[]
    for k1, k2 in itertools.combinations(v, 2):
        my_list.append([k1,k2, wordanalysis.euc_dist(v[k1], v[k2])])
    euc_list[k]=my_list
    f = "../Data/euclidean_distance"+str(k)+'.csv'
    # writing to csv file  
    with open(f, 'w') as csvfile:  
    # creating a csv writer object  
        csvwriter = csv.writer(csvfile)      
    # writing the fields  
        csvwriter.writerow(header)  
    # writing the data rows  
        csvwriter.writerows(my_list) 


