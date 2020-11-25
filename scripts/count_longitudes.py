#!/usr/bin/python
#-*- coding: utf-8 -*-
from collections import defaultdict
import numpy as np
import pandas as pd

def count_lengths(file):
    """Crear una lista con las longitudes del número de células
       por fila (cell file) """
    #read as list input file
    my_list = file.readlines()
    # you may also want to remove whitespace characters like `\n` at the end of each line
    my_list = [x.strip() for x in my_list]
    #Convertir todos los caracteres en mayúscula
    my_list = [x.upper() for x in my_list]
    #Create list of lengths
    length_cells = [] 
    #En cada elemento de la lista obtener la lingitud.
    for x in my_list:
        length_cells.append(len(x))  
    return length_cells

# Open the file and evaluate safe every line in a list of P. bracteatus
f=open('../Data/Pedilanthus/P_bracteatus/845_edited_cells.txt')

length_cells = count_lengths(f)
#Make a data frame with pandas
df = pd.DataFrame(length_cells, columns =['Number of cells'])
df.insert(0, "Sample", np.repeat("845", np.shape(df)[0]), True)
df.insert(0, "Species", np.repeat("Euphorbia bracteata", np.shape(df)[0]), True)
df.to_csv('../Data/Pedilanthus/P_bracteatus/bracteatus_length_filecells.csv', index = False)

#############Open files of P. calcaratus
f=open('../Data/Pedilanthus/P_calcaratus/892_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df = pd.DataFrame(length_cells, columns =['Number of cells'])
df.insert(0, "Sample", np.repeat("892", np.shape(df)[0]), True)
df.insert(0, "Species", np.repeat("Euphorbia calcarata", np.shape(df)[0]), True)

f=open('../Data/Pedilanthus/P_calcaratus/896_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df_1 = pd.DataFrame(length_cells, columns =['Number of cells'])
df_1.insert(0, "Sample", np.repeat("896", np.shape(df_1)[0]), True)
df_1.insert(0, "Species", np.repeat("Euphorbia calcarata", np.shape(df_1)[0]), True)

df.append(df_1)

df.to_csv('../Data/Pedilanthus/P_calcaratus/calcaratus_length_filecells.csv', index = False)

##############Open files of P. coalcomanensis
f=open('../Data/Pedilanthus/P_coalcomanensis/883_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df = pd.DataFrame(length_cells, columns =['Number of cells'])
df.insert(0, "Sample", np.repeat("883", np.shape(df)[0]), True)
df.insert(0, "Species", np.repeat("Euphorbia coalcomanensis", np.shape(df)[0]), True)

df.to_csv('../Data/Pedilanthus/P_coalcomanensis/coalcomanensis_length_filecells.csv', index = False)

##############Open files of P. colligata connatus
f=open('../Data/Pedilanthus/P_colligata/867_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df = pd.DataFrame(length_cells, columns =['Number of cells'])
df.insert(0, "Sample", np.repeat("867", np.shape(df)[0]), True)
df.insert(0, "Species", np.repeat("Euphorbia colligata", np.shape(df)[0]), True)

df.to_csv('../Data/Pedilanthus/P_colligata/colligata_length_filecells.csv', index = False)


##############Open files of P. cymbiferus
f=open('../Data/Pedilanthus/P_cymbiferus/979_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df = pd.DataFrame(length_cells, columns =['Number of cells'])
df.insert(0, "Sample", np.repeat("979", np.shape(df)[0]), True)
df.insert(0, "Species", np.repeat("Euphorbia cymbifera", np.shape(df)[0]), True)

df.to_csv('../Data/Pedilanthus/P_cymbiferus/cymbiferus_length_filecells.csv', index = False)


##############Open files of P. diazluna
f=open('../Data/Pedilanthus/P_diazluna/EPM10_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df = pd.DataFrame(length_cells, columns =['Number of cells'])
df.insert(0, "Sample", np.repeat("EPM10", np.shape(df)[0]), True)
df.insert(0, "Species", np.repeat("Euphorbia diazlunana", np.shape(df)[0]), True)

f=open('../Data/Pedilanthus/P_diazluna/EPM11_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df_1 = pd.DataFrame(length_cells, columns =['Number of cells'])
df_1.insert(0, "Sample", np.repeat("EPM11", np.shape(df_1)[0]), True)
df_1.insert(0, "Species", np.repeat("Euphorbia diazlunana", np.shape(df_1)[0]), True)

df.append(df_1)

f=open('../Data/Pedilanthus/P_diazluna/EPM12_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df_1 = pd.DataFrame(length_cells, columns =['Number of cells'])
df_1.insert(0, "Sample", np.repeat("EPM12", np.shape(df_1)[0]), True)
df_1.insert(0, "Species", np.repeat("Euphorbia diazlunana", np.shape(df_1)[0]), True)

df.append(df_1)

df.to_csv('../Data/Pedilanthus/P_diazluna/diazluna_length_filecells.csv', index = False)

#########################Files of E. finkii
f=open('../Data/Pedilanthus/P_finkii/917_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df = pd.DataFrame(length_cells, columns =['Number of cells'])
df.insert(0, "Sample", np.repeat("917", np.shape(df)[0]), True)
df.insert(0, "Species", np.repeat("Euphorbia finkii", np.shape(df)[0]), True)

df.to_csv('../Data/Pedilanthus/P_finkii/finkii_length_filecells.csv', index = False)

#########################Files of E. lomelii
f=open('../Data/Pedilanthus/P_macrocarpus/853_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df = pd.DataFrame(length_cells, columns =['Number of cells'])
df.insert(0, "Sample", np.repeat("853", np.shape(df)[0]), True)
df.insert(0, "Species", np.repeat("Euphorbia lomelii", np.shape(df)[0]), True)

df.to_csv('../Data/Pedilanthus/P_macrocarpus/macrocarpus_length_filecells.csv', index = False)


#########################Files of E. personata
f=open('../Data/Pedilanthus/P_personata/EPM7_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df = pd.DataFrame(length_cells, columns =['Number of cells'])
df.insert(0, "Sample", np.repeat("EPM7", np.shape(df)[0]), True)
df.insert(0, "Species", np.repeat("Euphorbia personata", np.shape(df)[0]), True)
#
f=open('../Data/Pedilanthus/P_personata/EPM9_edited_cells.txt')
length_cells = count_lengths(f)

df_1 = pd.DataFrame(length_cells, columns =['Number of cells'])
df_1.insert(0, "Sample", np.repeat("EPM9", np.shape(df_1)[0]), True)
df_1.insert(0, "Species", np.repeat("Euphorbia lomelii", np.shape(df_1)[0]), True)


df.append(df_1)
df.to_csv('../Data/Pedilanthus/P_personata/personata_length_filecells.csv', index = False)

#########################Files of E. peritropoides
f=open('../Data/Pedilanthus/P_peritropoides/974_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df = pd.DataFrame(length_cells, columns =['Number of cells'])
df.insert(0, "Sample", np.repeat("974", np.shape(df)[0]), True)
df.insert(0, "Species", np.repeat("Euphorbia peritropoides", np.shape(df)[0]), True)
#
df.to_csv('../Data/Pedilanthus/P_peritropoides/peritropoides_length_filecells.csv', index = False)

#########################Files of E. conzattii
f=open('../Data/Pedilanthus/P_pulchellus/971a_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df = pd.DataFrame(length_cells, columns =['Number of cells'])
df.insert(0, "Sample", np.repeat("971a", np.shape(df)[0]), True)
df.insert(0, "Species", np.repeat("Euphorbia conzattii", np.shape(df)[0]), True)
df.to_csv('../Data/Pedilanthus/P_pulchellus/pulchellus_length_filecells.csv', index = False)



#########################Files of E. tehuacana
f=open('../Data/Pedilanthus/P_tehuacanus/981_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df = pd.DataFrame(length_cells, columns =['Number of cells'])
df.insert(0, "Sample", np.repeat("981", np.shape(df)[0]), True)
df.insert(0, "Species", np.repeat("Euphorbia tehuacana", np.shape(df)[0]), True)
df.to_csv('../Data/Pedilanthus/P_tehuacanus/tehuacanus_length_filecells.csv', index = False)


#########################Files of E. tithymaloides
f=open('../Data/Pedilanthus/P_tithymaloides/EPM6_S2-1_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df = pd.DataFrame(length_cells, columns =['Number of cells'])
df.insert(0, "Sample", np.repeat("EPM6", np.shape(df)[0]), True)
df.insert(0, "Species", np.repeat("Euphorbia tithymaloides", np.shape(df)[0]), True)

f=open('../Data/Pedilanthus/P_tithymaloides/EPM5_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df_1 = pd.DataFrame(length_cells, columns =['Number of cells'])
df_1.insert(0, "Sample", np.repeat("EPM5", np.shape(df_1)[0]), True)
df_1.insert(0, "Species", np.repeat("Euphorbia tithymaloides", np.shape(df_1)[0]), True)


df.append(df_1)

df.to_csv('../Data/Pedilanthus/P_tithymaloides/tithymaloides_length_filecells.csv', index = False)

#########################Files of E. cyri
f=open('../Data/Pedilanthus/P_tomentellus/973_edited_cells.txt')
length_cells = count_lengths(f)
#Make a data frame with pandas
df = pd.DataFrame(length_cells, columns =['Number of cells'])
df.insert(0, "Sample", np.repeat("973", np.shape(df)[0]), True)
df.insert(0, "Species", np.repeat("Euphorbia cyri", np.shape(df)[0]), True)

df.to_csv('../Data/Pedilanthus/P_tomentellus/tomentellus_length_filecells.csv', index = False)



