#!/usr/bin/python
#-*- coding: utf-8 -*-

# Open the file and evaluate safe every line in a list
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
import numpy as np

def conduc_index(string): 
    """Function to transform a sequence of cells into
    a binary sequence of 0 and 1 based on conducitve properties"""
    c = ''
    string += ' '
    for z in range(0,len(string)):
        if string[z] == 'V':
            c += '1'    
        elif string[z] == 'F' and string[z+1] == 'V':    
            c += '1'
        elif string[z] == 'F' and string[z-1] == 'V':    
            c += '1'    
        elif string[z] == 'P' and string[z+1] == 'V':
            c += '1'
        elif string[z] == 'P' and string[z-1] == 'V':
            c += '1'
        elif string[z] == 'P' and string[z+1] == 'F':
            c += '0'
        elif string[z] == 'P' and string[z-1] == 'F':
            c += '0'
    return(c)

def storage_index(string):
    """Function to transform a sequence of cells into
    a binary sequence of 0 and 1 based on storage properties"""
    c = ''
    string += ' '
    for z in range(0,len(string)):
        if string[z] == 'P':    
            c += '1'
        elif string[z] == 'V':
            c += '0'    
        elif string[z] == 'F' and string[z+1] == 'P':
            c += '1'
        elif string[z] == 'F' and string[z-1] == 'P':
            c += '1'
        elif string[z] == 'F' and string[z+1] == 'V':
            c += '0'
        elif string[z] == 'F' and string[z-1] == 'V':
            c += '0'
    return(c)

def support_index(string):
    """Function to transform a sequence of cells into
    a binary sequence of 0 and 1 based on support properties"""
    c = ''
    string += ' '
    for z in range(0,len(string)):
        if string[z] == 'V':
            c += '0'
        elif string[z] == 'P':
            c += '0'
        elif string[z] == 'F' and string[z-1] == 'V':    
            c += '0'
        elif string[z] == 'F' and string[z+1] == 'V':    
            c += '0'
        elif string[z] == 'F':    
            c += '1'     
    return(c)

#Calcular el índice de homogeneidad
def homogenity_index(string):
    """"After recieving a binary sequence, calc the distribution 
    of 00, 01, 10, and 11"""
    n01 = 0
    n10 = 0
    n11 = 0
    n00 = 0
    n0 = 0
    n1 = 0
    for i in range(0,len(string)):
        if string[i:i+2] == '01':
            n01 += 1
        elif string[i:i+2] == '10':
            n10 += 1
        elif string[i:i+2] == '00':
            n00 += 1
        elif string[i:i+2] == '11':
            n11 += 1
    for i in range(0,len(string)): 
        if string[i]=='1':
            n1 += 1
        if string[i]=='0':
            n0 += 1
    if 0 in {n0,n1}:
        d=1
    elif 0 not in {n0,n1}:
        d=((n00*n11)-(n10*n01))/(n0*n1)    
    return(d)

##
def percentage01(string):  
    """FUnction to return the % of 0 and 1s in a sequence"""
    if string.count('0') == 0:
        return(1)
    elif string.count('1') == 0:
        return(0)
    else:
        ones = (string.count('1')/len(string))
        return(ones)

#First define the list of files that have the cell files
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
            
#Now dictoniary is ready for analysis
#Remove ray cells without changing the order of lineages
for values in filenames.values(): 
    values.append([]) 
    for i in range(0,len(values[1])):
        x = values[1][i].replace('R','')
        if x == '':
            pass
        else:
            values[2].append(x)

##Create an empty dictionary to save the binary files
# After creating dictionaries with binary sequence for each 
# index, calc the values of homogeneity
conductivity_code = {} #Empty dictionary to add conductivity values into
for i in range(0,len(filenames.keys())):
    conductivity_code[list(filenames.keys())[i]]=[] #add key element to dict
for k, v in filenames.items():
    for z in v[2]:
        conductivity_code[k].append(conduc_index(z)) 
#Make it for the storage index
storage_code = {} #Empty dictionary to add values into
for i in range(0,len(filenames.keys())):
    storage_code[list(filenames.keys())[i]]=[] #add key element to dict
for k, v in filenames.items():
    prueba_list=[]
    for z in v[2]:
        storage_code[k].append(storage_index(z))
#Make it for the SUPPORT index
support_code = {} #Empty dictionary to add values into
for i in range(0,len(filenames.keys())):
    support_code[list(filenames.keys())[i]]=[] #add key element to dict
for k, v in filenames.items():
    prueba_list=[]
    for z in v[2]:
        support_code[k].append(support_index(z))
conductivity_index = {} #Empty dictionary to add values into
for i in range(0,len(filenames.keys())):
    conductivity_index[list(filenames.keys())[i]]=[] 
for k, v in conductivity_code.items():
    prueba_list=[]
    for z in v:
         conductivity_index[k].append(homogenity_index(z)) 
storage_index = {} #Empty dictionary to add values into
for i in range(0,len(filenames.keys())):
    storage_index[list(filenames.keys())[i]]=[] 
for k, v in storage_code.items():
    prueba_list=[]
    for z in v:
         storage_index[k].append(homogenity_index(z)) 
#
support_index = {} #Empty dictionary to add values into
for i in range(0,len(filenames.keys())):
    support_index[list(filenames.keys())[i]]=[] 
for k, v in support_code.items():
    prueba_list=[]
    for z in v:
         support_index[k].append(homogenity_index(z)) 
#
my_lens = {} #empty dictionary to add length of cell files
for keys, values in filenames.items():
    with open(path+values[0]) as f:
        x = f.readlines()
        filenames.setdefault(keys,[]).append(x)
    values[2] = [x.strip() for x in values[2]]  
    values[2] = [x.upper() for x in values[2]]  
    my_lens[values[0]]=[]
    for s in values[2]:
        my_lens[values[0]].append(len(s))
#Create dict with values of percentage of 1's
conductivity_freq = {} #Empty dictionary to add values into
for i in range(0,len(filenames.keys())):
    conductivity_freq[list(filenames.keys())[i]]=[] 
for k, v in conductivity_code.items():
    prueba_list=[]
    for z in v:
         conductivity_freq[k].append(percentage01(z)) 
#For support
support_freq = {} #Empty dictionary to add values into
for i in range(0,len(filenames.keys())):
    support_freq[list(filenames.keys())[i]]=[] 
for k, v in support_code.items():
    prueba_list=[]
    for z in v:
         support_freq[k].append(percentage01(z)) 
                  
#####Create a pandas data frame
df = pd.DataFrame(columns=['sp'])
for k,v in storage_index.items():
    x = np.repeat(k,len(v))
    y = np.array(v)
    df = df.append(x.tolist())

indexall=[]
for k,v in storage_index.items():    
    for i in range(0,len(v)):
        indexall.append(v[i])
df["storage"]= indexall
indexall=[]
for k,v in conductivity_index.items():    
    for i in range(0,len(v)):
        indexall.append(v[i])
df["conductivity"]= indexall

indexall=[]
for k,v in support_index.items():    
    for i in range(0,len(v)):
        indexall.append(v[i])
df["support"]= indexall

indexall=[]
for k,v in my_lens.items():    
    for i in range(0,len(v)):
        indexall.append(v[i])
df["length"]= indexall

indexall=[]
for k,v in conductivity_freq.items():    
    for i in range(0,len(v)):
        indexall.append(v[i])
df["conductivityfreq"]= indexall    
indexall=[]

for k,v in support_freq.items():    
    for i in range(0,len(v)):
        indexall.append(v[i])
df["supportfreq"]= indexall    

indexall=[]
for k, v in filenames.items():
    for z in v[2]:
        indexall.append(z)
df["sequences"]= indexall       
        
df.to_csv('../Data/homogentiy_index.csv')

