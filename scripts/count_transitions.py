#!/usr/bin/python
#-*- coding: utf-8 -*-
from collections import defaultdict
# Asignar variable con valores de 0 a cada transición.
ff=0
fv=0
fp=0
vv=0
vp=0
vf=0
pp=0
pv=0
pf=0

# Open the fil and evaluate safe every line in a list
f=open(input('Ingresa el nombre del archivo ce células:'))
#Crear una lista con cada una de las líneas.
my_list = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
my_list = [x.strip() for x in my_list]

#Convertir todos los caracteres en mayúscula
my_list = [x.upper() for x in my_list]
# Loop 
z=0
y=0
for x in my_list: 
    z=0
    for z in range(0,len(my_list[y])):
        if my_list[y][z:z+2] == 'FF':
            ff+=1
        elif my_list[y][z:z+2] == 'FV': 
            fv+=1
        elif my_list[y][z:z+2] == 'FP':
            fp+=1
        elif my_list[y][z:z+2] == 'VV':
            vv+=1
        elif my_list[y][z:z+2] == 'VP':
            vp+=1
        elif my_list[y][z:z+2] == 'VF':
            vf+=1
        elif my_list[y][z:z+2] == 'PP':
            pp+=1
        elif my_list[y][z:z+2] == 'PV':
            pv+=1
        elif my_list[y][z:z+2] == 'PF':
            pf+=1
        z+=1
    y+=1

# Guardar en dos tablas los conteos de transiciones y además los nombres
transitions=[ff,fv,fp,vv,vp,vf,pp,pv,pf]
names=["F-F","F-V","F-P","V-V","V-P","V-F","P-P","P-V","P-F"]
#Hacer un loop para sacar las frecuencias.
frequencies=[]
#Ya está lista una lista vacía para hacer el loop y agregar lo que haga falta.
for i in transitions:
    x = i/sum(transitions[:])
    frequencies.append((i,x))
print(sum(transitions[:]))
#Imprimir las frecuecnias junto con sus nombres
for i in range(0,len(frequencies)):
    print(names[i], frequencies[i])

print(frequencies)

