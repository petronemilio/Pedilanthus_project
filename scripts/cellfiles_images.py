import pandas as pd
import numpy as np
from scipy import spatial
import scipy
import csv 

#Open concatenated file 
celltypeall = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/prueba_lineage.csv")
#
celltypeall['name'] = celltypeall['name'].apply(str)
###
lineage = list(set(celltypeall['Lineage']))
dictLists = dict((key, []) for key in lineage) #Create dictionary from each cell file
#num_cells = []
#
for keys in dictLists:
    # extract the rows
    subset = celltypeall[celltypeall['Lineage'] == keys]
    dictLists[keys].append(''.join(subset['name'].tolist())) #Extract the cell types from each cell file
    
def replace_all(text, dic):
    ###This function replace characters 
    ###in strings ###
    for i, j in dic.items():
        text = text.replace(i, j)
    return text
###Define dictionary with the characters to replace
char_to_replace = {'1': 'F','2': 'V',
                   '3': 'P','4': 'R'}
#Create new dictionary to replace numbers with letters of the different cell types
new_dict ={}
for key, values in dictLists.items():
    print(key,replace_all(values[0], char_to_replace))
    new_dict.update({key: replace_all(values[0], char_to_replace)})
    #key.append(replace_all(values[0], char_to_replace))

#Save data frame 
file = open("../Data/Images/P_tithymaloides/EPM6/lineage_transformed.txt"
            ,"w")
 
for value in new_dict.values(): 
    file.write('{}\n'.format(value))
    
file.close()
