#!/usr/bin/python
#-*- coding: utf-8 -*-

# Open the file and evaluate safe every line in a list
f=open(input('Ingresa el nombre del archivo ce células:'))
#Crear una lista con cada una de las líneas.
my_list = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
my_list = [x.strip() for x in my_list]

#Convertir todos los caracteres en mayúscula
my_list = [x.upper() for x in my_list]

#Crear las listas donde se guardaran las secuencias en código binario
#vivo_binarias para las células vivas y muertas
vivo_binarias=[] #P=1, V y F= 0
#cond_binarias=[] Conductor vs no conductor
cond_binarias=[] #V=1, P y F= 0
#sosten 
sost_binarias=[] #F=1, P y V= 0

#Con un loop se trasnforman las cadenas a binarias dependiendo de el 
#indice 
y=0
for i in my_list:
    z=0
    x='' # indice de vivas/muertas
    a='' # indice de conductoras/noconductoras
    s='' # indice de sosten/nososten	
    for z in range(0,len(my_list[y])):
        if my_list[y][z] == 'F':    
            x += '0'
            a += '0'
            s += '1'
        elif my_list[y][z] == 'V':
            x += '0'
            a += '1'
            s += '0'
        elif my_list[y][z] == 'P':   
            x += '1'
            a += '0'
            s += '0'
        z += 1
    vivo_binarias.append(x)
    cond_binarias.append(a)
    sost_binarias.append(s)
    y+=1
#Calcular el índice de homogeneidad
def homogenity_index(string):
    n01=0
    n11=0
    n0=0
    n1=1
    for i in range(0,len(string)):
        if string[i:i+2] == '01':
            n01+=1
        elif string[i:i+2] == '10':
            n01+=1
        elif string[i:i+2] == '00':
            n11+=1
        elif string[i:i+2] == '11':
            n11+=1
    for i in range(0,len(string)): 
        if string[i]=='1':
            n0+=1
        if string[i]=='0':
            n1+=1
    if 0 in {n0,n1}:
        d=1
    elif 0 not in {n0,n1}:
        d=(n01-n11)/(n0*n1)
    return(d)

vivo_values=[]
cond_values=[]
sost_values=[]
for f in vivo_binarias:
    vivo_values.append(homogenity_index(f))
for i in cond_binarias:
    cond_values.append(homogenity_index(i))
for i in sost_binarias:
    sost_values.append(homogenity_index(s))


fmt = '{:<8}{:<30}{}'

print(fmt.format('', 'Living-Death', 'Conductive-Non conductive'))
for i, (name, grade) in enumerate(zip(vivo_values, cond_values)):
    print(fmt.format(i, name, grade))
