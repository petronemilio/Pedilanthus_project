import string
import random

list_strings=list(string.ascii_lowercase)
list_strings.append(" ")
x = "methinks it is like a weasel"
i = ""
z=0
while i != x: 
    c = random.choice(list_strings) 
    if c == x[z]:
        i = i + c
        z +=1
            
