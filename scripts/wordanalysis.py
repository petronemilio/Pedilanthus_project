#!/usr/bin/python
#-*- coding: utf-8 -*-
###Here we have functions to count words and analyse cell files from secondary xylem derived cells.
import math
def count_cells(file):
    """
    Returns list with number of cells per wood sample.
    """
    #read as list input file
    my_list = file.readlines()
    file = [x.strip() for x in my_list]
    vasos = 0
    fibras = 0
    parenquima = 0
    radio = 0
    cell_counts={}    
    for x in my_list: 
        vasos += x.count('V')
        fibras += x.count('F')
        parenquima += x.count('P')
        radio += x.count('R')
    cell_counts["Vasos"]= vasos
    cell_counts["Fibras"]= fibras
    cell_counts["Parénquima"]= parenquima
    cell_counts["Radios"]= radio
    return cell_counts

def count_lengths(file):
    """
    Returns a list with the length of each cell file of a 
    sample wood
    """
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

#Get words regresa un diccionario con las palabras y los conteos.
#Se puede aprovechar para usarlo con multiples conteos
def get_words(file, x):
    """Determine the number of different words in 
    a given string. It needs two arguments: 
    1)a given string
    2) the length of the k-mer to determine"""
    #Create empty set to add the different words
    words = {}
    #go trhough the string:
    for z in range(0,len(file)):
        word = file[z:z+x]
        if len(word) == x:
            words[word] = words.get(word,0)
            words[word] = words[word] +1
        else:
            pass
    return(words)

def get_all_words(array, n):
    """Determine the number of different words in 
    a given array. It uses the get_words function
    It also need a k-number for the length of words"""
    #Create empty set to add the different words
    words = {}
    #iterate trhough the files:
    for i in array:
        x = get_words(i, n) 
        for keys,values in x.items():
            words[keys] = words.get(keys,0)
            words[keys] = words[keys] + values 
        #aquí debería de venir una forma de agregar cada 
        #diccionario de x a words
    return(words)
#Se cambió la función..... en lugar..
#Parece que las dos funciones funcionan.#Utilizarlas con los datos

def euc_dist(x, y):
    """Function to calculate the euclidian distance between
    two different dictionaries with word counts."""
    ss = 0 #sum of squares between words
    for key in x:
        if key in y:
            ss += (x[key] - y[key])**2
        else:
            ss += (x[key])**2
    for key in y:
        if key not in x:
            ss += (y[key])**2
    distance = math.sqrt(ss)
    
    return distance
#Como sumar las palabras que no aparecen                     

