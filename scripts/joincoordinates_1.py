import pandas as pd
import numpy as np
from scipy import spatial
import scipy
import csv 

cellfile1 = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión1/results_celltype_both.csv")
cellfile2 = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión2/results_celltype_both.csv")
cellfile4 = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión4/results_celltype_both.csv")
cellfile6 = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión6/results_celltype_both.csv")
cellfile6_2 = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/results_celltype_both_2.csv")
cellfile7 = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión7/results_celltype_both.csv")
cellfile9 = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión9/results_celltype_both.csv")
cellfile11 = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión11/results_celltype_both.csv")
cellfile12 = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión12/results_celltype_both.csv")
cellfile14 = pd.read_csv("../Data/Images/P_tithymaloides/EPM6/10x/fusión14/results_celltype_both.csv")


cellfile = pd.concat([cellfile1,cellfile2,cellfile4,cellfile6,cellfile6_2,cellfile7,
           cellfile9,cellfile11,cellfile12,cellfile14])
           
#cellfile.insert(0, 'New_ID', range(880, 880 + len(df)))
newindex = []
z = 1
#trying to add a cell identifier to the concats in order to have the ordered sequence of cell lineages. 
cellfile.to_csv('../../../media/emiliopetronem/ADATA HD710 PRO/Pedilanthus/E_tithymaloides/EPM6/prueba.csv')
#File edited online

