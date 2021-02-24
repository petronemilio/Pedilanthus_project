#!/usr/bin/python
# -*- coding: utf-8 -*-
##Hay que decirle a python que utilice el coding utf-8 y no el ASCII character
import os
import shutil
import tempfile
import canvasvg
import re
import turtle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def nextpoint(actualpoint, character):
    """Recieves a list with x,y positions and a character 
    with info of the next coordinates"""
    coordenadas ={ 'V':(-250,250), 
                  'F':(-250,-250),'P':(250,250),'R':(250,-250)} 
    nextpoint= actualpoint
    nextpoint[0] = (coordenadas[character][0] + actualpoint[0])/2
    nextpoint[1] = (coordenadas[character][1] + actualpoint[1])/2
    return(nextpoint)

def chaosmap(fila):
    steps = np.zeros((len(fila),2))
    x = fila
    for i in range(1,len(x)):
        steps[i] = nextpoint(steps[i-1], x[i-1])
    return(steps)

#cells845 = []
#cells883 = []
#lsys = []
#rlsys = []
#with open('../Data/Cell_files_data/883_edited_cells.txt') as f:
 #   for x in f:
  #      cells883.append(x)
#cells883 = [x.strip() for x in cells883]
#cells883 = [x.upper() for x in cells883]
#cells883iterator = []

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

chaos_representation = {} #Empty dictionary to add values into
for values in filenames.values():
    chaos_representation[values[0]]=[]
    for z in range(len(values[1])):
        chaos_representation[values[0]].append(chaosmap(values[1][z])) 

cols = ['X','Y', 'count','Sample']
dat = pd.DataFrame(columns = cols)
for k,v in chaos_representation.items():
    for y in v: #iterate through cell files for each sample
        x = 0   #create a counter to iterate trhough each cell in cell files  
        for z in range(x,len(y)):
            dat = dat.append({'X': y[z][0],'Y': y[z][1],'count':x, 'Sample': k}, ignore_index=True)
            x += 1
dat.to_csv('../Data/chaos_representation.csv')

#figures = {} #Empty dictionary to add values into
#for keys, values in chaos_representation.items():
 #   figures[keys]=[]
  #  for z in range(0,len(values)):
   #     x, y = values[z].T
    #    figures[keys]=plt.scatter(x,y)
     #   print(len(x),len(y)) 
#for keys, values in figures.items():
 #   values.figure.savefig('../Figures/'+ str(keys)+'.pdf')


####### trying turtle

#height = 360
#width = 360

#wn = turtle.Screen()
#wn.screensize(width, height)
#wn.bgcolor("white")

#ojoche = turtle.Turtle()
#ojoche.shape('circle')
#ojoche.turtlesize(0.3)
#ojoche.penup()
#ojoche.pensize(2)  
#ojoche.speed('fastest')

#for i in rlsysiterator:
#    ojoche.goto(0,0)
 #   for z in i:
  #      ojoche.goto( turtle.pos() + (z[0],z[1]) )
   #     ojoche.pendown()
    #    ojoche.stamp()
     #   ojoche.penup()
      #  ojoche.goto(0,0)
    

#wn.exitonclick()                # wait for a user click on the canvas

 
