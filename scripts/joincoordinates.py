import pandas as pd
import numpy as np
from scipy import spatial
import scipy
import csv 
###Temporal loading of each file
celltype = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión1/results_celltype.csv")
cellfile = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión1/results_file.csv")
#####
pairs1 = list(tuple(zip(celltype.X, celltype.Y)))
pairs2 = list(tuple(zip(cellfile.X, cellfile.Y)))
min_distances = []
closest_pairs = []
names = []
for i in pairs2:
    min_dist = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').min()
    index_min = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').argmin()
    min_distances.append(min_dist)
    closest_pairs.append(celltype.loc[index_min, ['X', 'Y']])
    names.append(celltype.loc[index_min, 'Type'])

cellfile['min_distance'] = min_distances
cellfile['closest_pairs'] = [tuple(i.values) for i in closest_pairs]
cellfile['name'] = names
######
cellfile.to_csv('../Data/Images/P_tithymaloides/EPM6/10x/fusión1/results_celltype_both.csv', index=False) 
######write second file
celltype = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión2/results_celltype.csv")
cellfile = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión2/results_file.csv")
#####
pairs1 = list(tuple(zip(celltype.X, celltype.Y)))
pairs2 = list(tuple(zip(cellfile.X, cellfile.Y)))
min_distances = []
closest_pairs = []
names = []
for i in pairs2:
    min_dist = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').min()
    index_min = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').argmin()
    min_distances.append(min_dist)
    closest_pairs.append(celltype.loc[index_min, ['X', 'Y']])
    names.append(celltype.loc[index_min, 'Type'])

cellfile['min_distance'] = min_distances
cellfile['closest_pairs'] = [tuple(i.values) for i in closest_pairs]
cellfile['name'] = names
######
cellfile.to_csv('../Data/Images/P_tithymaloides/EPM6/10x/fusión2/results_celltype_both.csv', index=False) 

######OPEN AND WRITE OTHER FUSION
celltype = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión4/results_celltype.csv")
cellfile = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión4/results_file.csv")
#####
pairs1 = list(tuple(zip(celltype.X, celltype.Y)))
pairs2 = list(tuple(zip(cellfile.X, cellfile.Y)))
min_distances = []
closest_pairs = []
names = []
for i in pairs2:
    min_dist = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').min()
    index_min = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').argmin()
    min_distances.append(min_dist)
    closest_pairs.append(celltype.loc[index_min, ['X', 'Y']])
    names.append(celltype.loc[index_min, 'Type'])

cellfile['min_distance'] = min_distances
cellfile['closest_pairs'] = [tuple(i.values) for i in closest_pairs]
cellfile['name'] = names
######
cellfile.to_csv('../Data/Images/P_tithymaloides/EPM6/10x/fusión4/results_celltype_both.csv', index=False)

######OPEN ANOTHER SET OF FILES 6
celltype = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión6/results_celltype.csv")
cellfile = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión6/results_file.csv")
#####
pairs1 = list(tuple(zip(celltype.X, celltype.Y)))
pairs2 = list(tuple(zip(cellfile.X, cellfile.Y)))
min_distances = []
closest_pairs = []
names = []
for i in pairs2:
    min_dist = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').min()
    index_min = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').argmin()
    min_distances.append(min_dist)
    closest_pairs.append(celltype.loc[index_min, ['X', 'Y']])
    names.append(celltype.loc[index_min, 'Type'])

cellfile['min_distance'] = min_distances
cellfile['closest_pairs'] = [tuple(i.values) for i in closest_pairs]
cellfile['name'] = names
######
cellfile.to_csv('../Data/Images/P_tithymaloides/EPM6/10x/fusión6/results_celltype_both.csv', index=False)
######OPEN ANOTHER SET OF FILES 6-2
celltype = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión6/results_celltype_2.csv")
cellfile = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión6/results_file_2.csv")
#####
pairs1 = list(tuple(zip(celltype.X, celltype.Y)))
pairs2 = list(tuple(zip(cellfile.X, cellfile.Y)))
min_distances = []
closest_pairs = []
names = []
for i in pairs2:
    min_dist = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').min()
    index_min = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').argmin()
    min_distances.append(min_dist)
    closest_pairs.append(celltype.loc[index_min, ['X', 'Y']])
    names.append(celltype.loc[index_min, 'Type'])

cellfile['min_distance'] = min_distances
cellfile['closest_pairs'] = [tuple(i.values) for i in closest_pairs]
cellfile['name'] = names
######
cellfile.to_csv('EPM6/10x/fusión6/results_celltype_both_2.csv', index=False)
######OPEN ANOTHER SET OF FILES 7
celltype = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión7/results_celltype.csv")
cellfile = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión7/results_file.csv")
#####
pairs1 = list(tuple(zip(celltype.X, celltype.Y)))
pairs2 = list(tuple(zip(cellfile.X, cellfile.Y)))
min_distances = []
closest_pairs = []
names = []
for i in pairs2:
    min_dist = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').min()
    index_min = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').argmin()
    min_distances.append(min_dist)
    closest_pairs.append(celltype.loc[index_min, ['X', 'Y']])
    names.append(celltype.loc[index_min, 'Type'])

cellfile['min_distance'] = min_distances
cellfile['closest_pairs'] = [tuple(i.values) for i in closest_pairs]
cellfile['name'] = names
######
cellfile.to_csv('../Data/Images/P_tithymaloides/EPM6/10x/fusión7/results_celltype_both.csv', index=False)
######OPEN ANOTHER SET OF FILES 9
celltype = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión9/results_celltype.csv")
cellfile = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión9/results_file.csv")
#####
pairs1 = list(tuple(zip(celltype.X, celltype.Y)))
pairs2 = list(tuple(zip(cellfile.X, cellfile.Y)))
min_distances = []
closest_pairs = []
names = []
for i in pairs2:
    min_dist = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').min()
    index_min = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').argmin()
    min_distances.append(min_dist)
    closest_pairs.append(celltype.loc[index_min, ['X', 'Y']])
    names.append(celltype.loc[index_min, 'Type'])

cellfile['min_distance'] = min_distances
cellfile['closest_pairs'] = [tuple(i.values) for i in closest_pairs]
cellfile['name'] = names
######
cellfile.to_csv('../Data/Images/P_tithymaloides/EPM6/10x/fusión9/results_celltype_both.csv', index=False)
######OPEN ANOTHER SET OF FILES 11
celltype = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión11/results_celltype.csv")
cellfile = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión11/results_file.csv")
#####
pairs1 = list(tuple(zip(celltype.X, celltype.Y)))
pairs2 = list(tuple(zip(cellfile.X, cellfile.Y)))
min_distances = []
closest_pairs = []
names = []
for i in pairs2:
    min_dist = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').min()
    index_min = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').argmin()
    min_distances.append(min_dist)
    closest_pairs.append(celltype.loc[index_min, ['X', 'Y']])
    names.append(celltype.loc[index_min, 'Type'])

cellfile['min_distance'] = min_distances
cellfile['closest_pairs'] = [tuple(i.values) for i in closest_pairs]
cellfile['name'] = names
######
cellfile.to_csv('../Data/Images/P_tithymaloides/EPM6/10x/fusión11/results_celltype_both.csv', index=False)
#
######OPEN ANOTHER SET OF FILES 12
celltype = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión12/results_celltype.csv")
cellfile = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión12/results_file.csv")
#####
pairs1 = list(tuple(zip(celltype.X, celltype.Y)))
pairs2 = list(tuple(zip(cellfile.X, cellfile.Y)))
min_distances = []
closest_pairs = []
names = []
for i in pairs2:
    min_dist = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').min()
    index_min = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').argmin()
    min_distances.append(min_dist)
    closest_pairs.append(celltype.loc[index_min, ['X', 'Y']])
    names.append(celltype.loc[index_min, 'Type'])

cellfile['min_distance'] = min_distances
cellfile['closest_pairs'] = [tuple(i.values) for i in closest_pairs]
cellfile['name'] = names
######
cellfile.to_csv('../Data/Images/P_tithymaloides/EPM6/10x/fusión12/results_celltype_both.csv', index=False)
######OPEN ANOTHER SET OF FILES 14
celltype = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión14/results_celltype.csv")
cellfile = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión14/results_file.csv")
#####
pairs1 = list(tuple(zip(celltype.X, celltype.Y)))
pairs2 = list(tuple(zip(cellfile.X, cellfile.Y)))
min_distances = []
closest_pairs = []
names = []
for i in pairs2:
    min_dist = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').min()
    index_min = scipy.spatial.distance.cdist([i], pairs1, metric='euclidean').argmin()
    min_distances.append(min_dist)
    closest_pairs.append(celltype.loc[index_min, ['X', 'Y']])
    names.append(celltype.loc[index_min, 'Type'])

cellfile['min_distance'] = min_distances
cellfile['closest_pairs'] = [tuple(i.values) for i in closest_pairs]
cellfile['name'] = names
######
cellfile.to_csv('../Data/Images/P_tithymaloides/EPM6/10x/fusión14/results_celltype_both.csv', index=False)
