#!/usr/bin/python
#-*- coding: utf-8 -*-
###Here we have functions to count words and analyse cell files from secondary xylem derived cells.
import math
import collections
 
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


def lzcomprdist(string):
    """This function returns the Lempl ziv
    compression distance of each file of cells"""
    keys_dict = {}
    ind = 0
    inc = 1
    while True:
        if not (len(string) >= ind+inc):
            break
        sub_str = string[ind:ind + inc]
    #print(sub_str,ind,inc)
        if sub_str in keys_dict:
            inc += 1
        else:
            keys_dict[sub_str] = 0
            ind += inc
            inc = 1
    return(len(keys_dict))                  


def estimate_shannon_entropy(file):
    """This function returns Shannon entropy value of a string.
    collections.Counter makes a dictionary with counts of each letter in string.
    """    
    m = len(file)
    cells = collections.Counter([tmp_cell for tmp_cell in file])  
 
    shannon_entropy_value = 0
    for cell in cells:
        # number of residues
        n_i = cells[cell]
        # n_i (# residues type i) / M (# residues in column)
        p_i = n_i / float(m)
        entropy_i = p_i * (math.log(p_i, 2))
        shannon_entropy_value += entropy_i
 
    return shannon_entropy_value * -1

def estimate_fusiform_to_ray_transitions(string):
    """This function count the number of times a ray to fusiform
    transition in the files """
    R = 0
    for i in range(0,len(string)):
        if string[i:i+2] == 'RF':
            R += 1
        elif string[i:i+2] == 'RV':
            R += 1
        elif string[i:i+2] == 'RP':
            R += 1
    return(R)

    
def estimate_shannon_entropy_window(file):
    m = len(file)
    shannon_entropy_window=[]
    for i in range(0,int(m/5)):
        x = file[(i*5):((i*5)+15)]
        m_w = len(x)
        if m_w < 10:
            break
        cells = collections.Counter([tmp_cell for tmp_cell in x])
        #    
        shannon_entropy_value = 0
        for cell in cells:
            # number of residues
            n_i = cells[cell]
            # n_i (# residues type i) / M (# residues in column)
            p_i = n_i / float(m_w)
            entropy_i = p_i * (math.log(p_i, 2))
            shannon_entropy_value += entropy_i
        shannon_final = shannon_entropy_value * -1
        shannon_entropy_window.append(shannon_final)
    return shannon_entropy_window
