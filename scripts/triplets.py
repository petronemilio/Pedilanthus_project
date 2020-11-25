#!/usr/bin/python
#-*- coding: utf-8 -*-
from collections import defaultdict
# Asignar variable con valores de 0 a cada transición.
fff=0
ffv=0
ffp=0
fvf=0
fvv=0
fvp=0
fpf=0
fpv=0
fpp=0
vff=0
vfv=0
vfp=0
vvf=0
vvv=0
vvp=0
vpf=0
vpv=0
vpp=0
pff=0
pfv=0
pfp=0
pvf=0
pvv=0
pvp=0
ppf=0
ppv=0
ppp=0

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
        if my_list[y][z:z+3] == 'FFF':
            fff+=1
        elif my_list[y][z:z+3] == 'FFV': 
            ffv+=1
        elif my_list[y][z:z+3] == 'FFP': 
            ffp+=1
        elif my_list[y][z:z+3] == 'FVF': 
            fvf+=1
        elif my_list[y][z:z+3] == 'FVV': 
            fvv+=1
        elif my_list[y][z:z+3] == 'FVP': 
            fvp+=1
        elif my_list[y][z:z+3] == 'FPF': 
            fpf+=1
        elif my_list[y][z:z+3] == 'FPV': 
            fpv+=1
        elif my_list[y][z:z+3] == 'FPP': 
            fpp+=1
        #
        elif my_list[y][z:z+3] == 'VFF':
            vff+=1
        elif my_list[y][z:z+3] == 'VFV': 
            vfv+=1
        elif my_list[y][z:z+3] == 'VFP': 
            vfp+=1
        elif my_list[y][z:z+3] == 'VVF': 
            vvf+=1
        elif my_list[y][z:z+3] == 'VVV': 
            vvv+=1
        elif my_list[y][z:z+3] == 'VVP': 
            vvp+=1
        elif my_list[y][z:z+3] == 'VPF': 
            vpf+=1
        elif my_list[y][z:z+3] == 'VPV': 
            vpv+=1
        elif my_list[y][z:z+3] == 'VPP': 
            vpp+=1
        #
        elif my_list[y][z:z+3] == 'PFF':
            pff+=1
        elif my_list[y][z:z+3] == 'PFV': 
            pfv+=1
        elif my_list[y][z:z+3] == 'PFP': 
            pfp+=1
        elif my_list[y][z:z+3] == 'PVF': 
            pvf+=1
        elif my_list[y][z:z+3] == 'PVV': 
            pvv+=1
        elif my_list[y][z:z+3] == 'PVP': 
            pvp+=1
        elif my_list[y][z:z+3] == 'PPF': 
            ppf+=1
        elif my_list[y][z:z+3] == 'PPV': 
            ppv+=1
        elif my_list[y][z:z+3] == 'PPP': 
            ppp+=1
        z+=1
    y+=1

# Guardar en dos tablas los conteos de transiciones y además los nombres
triplets=[fff,ffv,ffp,fvf,fvf,fvp,fpf,fpv,fpp,
             vff,vfv,vfp,vvf,vvf,vvp,vpf,vpv,vpp,
	     pff,pfv,pfp,pvf,pvf,pvp,ppf,ppv,ppp]

names=["FFF","FFV","FFP","FVF","FVV","FVP","FPF","FPV","FPP",
	"VFF","VFV","VFP","VVF","VVV","VVP","VPF","VPV","VPP",
	"PFF","PFV","PFP","PVF","PVV","PVP","PPF","PPV","PPP"]
#Hacer un loop para sacar las frecuencias.
frequencies=[]
#Ya está lista una lista vacía para hacer el loop y agregar lo que haga falta.
for i in triplets:
    x = i/sum(triplets[:])
    frequencies.append((i,x))
print(sum(triplets[:]))
#Imprimir las frecuecnias junto con sus nombres
for i in range(0,len(frequencies)):
    print(names[i], frequencies[i])

print(frequencies)

