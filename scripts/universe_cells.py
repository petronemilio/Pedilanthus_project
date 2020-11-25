#!/usr/bin/python
# -*- coding: utf-8 -*-
##Hay que decirle a python que utilice el coding utf-8 y no el ASCII character

potencia= input("Introduce la longitud de la cadena a la que deseas calcular su universo:")
#Hay que convertir el input de str a int
potencia_int= int(potencia)

cadenas = 2 ** potencia_int

### En la función combinaciones, para un número dado se generan las
# cadenas que se pueden construir de esa longitud. 
def combinaciones(n):
    for i in range(2**n):
        s = bin(i)[2:] #Convertir el entero en valor binario. El [2:]es para dejar solo el binario 
        s = "0" * (n-len(s)) + s #Sacar en orden las cadenas binarias
        yield s

#El archivo se va a escribir a un archivo llamado cadenas.txt
f= open("cadenas.txt","w+")
#En el siguiente for loop se hace que salgan las cadenas posibles menores
#a la potencia ingresada.
every_iteration=0
for i in range(potencia_int):
    print >>f, list(combinaciones(every_iteration))    
    every_iteration= every_iteration+1

#Hacer una rutina para calcular las cadenas de 0 a 20 y contar el número de unos
