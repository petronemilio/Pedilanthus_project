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
## Script to count total number of words at different k-mer lengths 
###Load all files using os and re
path = '../Data/Cell_files_data/ConvergeAssOtherLineage/'
files = os.listdir(path)
filenames = {} #Create a dictionary to save paths for all data cells
for i in files:
    m = re.search(r'.*[^_edited_cells_NotConverge.txt]',i)
    filenames[m.group()] = [i]

#Loop to load the cell lineages and append them in filenames within a dictionary
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
#Remove ray cells without changing the order of fusiform derived cells.
for values in filenames.values(): 
    values.append([])
    for i in range(0,len(values[1])):
        x = values[1][i].replace('R','')
        if x == '':
            pass
        else:
            values[2].append(x)

#################  Create a dictionary to save the word counting #########
word_diversity = {} #Empty dictionary to add values of word counts into
for i in range(1,35): #loop to add words at k-mer length from 1 to 35
    word_diversity[i]={} 
    for values in filenames.values():   
        word_diversity[i][values[0]] = wordanalysis.get_all_words(values[2], i)

###Write matrix output from df
val = 1
for x in word_diversity.values():
    wordcount = pd.DataFrame.from_dict(x, orient='index')
    wordcount.to_csv('../Data/word_counts_all/wordcounts'+ str(val)+'.csv', index = True)    
    val += 1
####Write the total number of words at each length for samples:
##Word counts/features
word_features = []
for y in word_diversity.values():
    for a,b in y.items():
        word_features.append([a,len(b)])
#Convert the list into dataframe
data = pd.DataFrame(word_features, columns=['file','NumberOfWords'])
data.to_csv('../Data/word_counts_all/wordcounts_all.csv')

####Repeat but counting only words that appear more than once: #############
#
word_diversity = {} #Empty dictionary to add values into
for i in range(1,35):
    word_diversity[i]={}
    for values in filenames.values(): 
        word_diversity[i][values[0]] = wordanalysis.get_all_words_more_than_one(values[2], i)

###Write matrix output from df
val = 1
for x in word_diversity.values():
    wordcount = pd.DataFrame.from_dict(x, orient='index')
    wordcount.to_csv('../Data/word_counts_morethanone/wordcounts'+ str(val)+'.csv', index = True)    
    val += 1
####Write the total number of words at each length for samples:
##Word counts/features
word_features = []
for y in word_diversity.values():
    for a,b in y.items():
        word_features.append([a,len(b)])
#Convert the list into dataframe
data = pd.DataFrame(word_features, columns=['file','NumberOfWords'])
data.to_csv('../Data/word_counts_morethanone/wordcounts_all.csv')

