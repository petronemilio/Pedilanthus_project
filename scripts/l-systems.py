import random
import matplotlib.pyplot as plt
import numpy as np
import turtle
import math

###Script to create different types of L-systems

def applyRules(lhch):
    """Apply stochastic rules using random distribution
    for values between 0 and 2 to apply rul1 or rule2 or
    rule3 """
    rhstr = ""
    s = np.random.randint(0,3,1)
    if lhch == 'A' and s == 0:
        rhstr = 'AF'   # Rule 1
    elif lhch == 'A' and s == 1:
        rhstr = 'AV' #Rule 2
    elif lhch == 'A' and s == 2:
        rhstr = 'AP' #Rule 3
    elif lhch == 'A':
        rhstr = 'B'  # Rule 2
    else:
        rhstr = lhch    # no rules apply so keep the character

    return rhstr

def processString(oldStr):
    newstr = ""
    for ch in oldStr:
        newstr = newstr + applyRules(ch)
    return newstr

def createLSystem(numIters,axiom):
    startString = axiom
    endString = ""
    for i in range(numIters):
        endString = processString(startString)
        startString = endString
    return endString

##Apply the previuos functions to get an array of files
#in which cells had equal probability of formation each cycle
lengths_array = np.random.randint(10,150,120)
random_files = []
for i in lengths_array:
    random_files.append(createLSystem(i,'A'))
#Remove A and B
for x, y in enumerate(random_files): 
    random_files[x] = y.replace('B', '') 
for x, y in enumerate(random_files):     
    random_files[x] = y.replace('A', '') 
    
####    
with open('../Data/Cell_files_data/probLsystem_edited_cells.txt', 'w') as output:
    for row in random_files:
        output.write(str(row) + '\n')


####L-system with context sensitivity for average trasnition of xeric habit
#F-F	F-V	    F-P	    0.9288	0.0017	0.0695
#V-V    V-P     V-F     0.2784	0.3658	0.3558
#P-P	P-V	    P-F     0.2247	0.1427	0.6326

rules_ctsens = { "A" : "AB" , "B" : ["V","F","P"]}
def iteratex(axioms, rules):
    new_text = []
    parsero = list(axioms)
    parsero.append('')
    for i, item in enumerate(parsero):
        if item == '':
            pass
        #If string is just AB it will generate V,F or P with equal probability
        elif item in rules and parsero[i+1]=='':
            s = np.random.randint(0,3,1)
            if s == 0:
                new_text.append(rules[item][0])
            elif s == 1:
                new_text.append(rules[item][1])
            else:
                new_text.append(rules[item][2])
        #Contexte sensitivity. If V is the last derivative produced 
        # probabilities are different
        elif item in rules and parsero[i+1]=='V':  # 0.2784	0.3658	0.3558
            s = np.random.uniform(low=0, high=1, size=None)
            if s < 0.2784:
                new_text.append(rules[item][0])
            elif 0.2784 <= s < 0.6442:
                new_text.append(rules[item][2])  
            else:    
                new_text.append(rules[item][1])     
        #Contexte sensitivity. If F is the last derivative produced 
        # probabilities are different
        elif item in rules and parsero[i+1]== 'F':
            s = np.random.uniform(low=0, high=1, size=None)
            if s < 0.9288:   #0.9288	0.0017	0.0695
                new_text.append(rules[item][1])
            elif 0.9288 <= s < 0.9305:
                new_text.append(rules[item][0])  
            else:    
                new_text.append(rules[item][2])
        #Contexte sensitivity. If P is the last derivative produced 
        # probabilities are different        
        elif item in rules and parsero[i+1]== 'P': #0.2247	0.1427	0.6326
            s = np.random.uniform(low=0, high=1, size=None)
            if s < 0.2247 :
                new_text.append(rules[item][2])
            elif 0.2247 <= s < 0.3674 :
                new_text.append(rules[item][0])  
            else:    
                new_text.append(rules[item][1])  
                
        elif item in rules:
            new_text.append(rules[item])
        
        else:
            new_text.append(item)
    return ''.join(new_text)

def clsystemx(axioms, rules, iteration):
    x = axioms
    for i in range(iteration):
        x=iteratex(x, rules)
    return ''.join(x)

##Apply the previuos functions to get an array of files
#in which cells had equal probability of formation each cycle
lengths_array = np.random.randint(10,150,120)
random_files = []
for i in lengths_array:
    random_files.append(clsystemx('AB',rules_ctsens,i))
#Remove A and B
for x, y in enumerate(random_files): 
    random_files[x] = y.replace('B', '') 
for x, y in enumerate(random_files):     
    random_files[x] = y.replace('A', '') 
    
with open('../Data/Cell_files_data/contextxericLsystem_edited_cells.txt', 'w') as output:
    for row in random_files:
        output.write(str(row) + '\n')

###########     Repeat for the mesic transitions 
####L-system with context sensitivity for average trasnition of xeric habit
#F-F	F-V	    F-P	    0.8151	0.0015	0.1834
#V-V    V-P     V-F     0.4804	0.3513	0.1683
#P-P	P-V	    P-F     0.0926	0.0418	0.8656

def iteratem(axioms, rules):
    new_text = []
    parsero = list(axioms)
    parsero.append('')
    for i, item in enumerate(parsero):
        if item == '':
            pass
        #If string is just AB it will generate V,F or P with equal probability
        elif item in rules and parsero[i+1]=='':
            s = np.random.randint(0,3,1)
            if s == 0:
                new_text.append(rules[item][0])
            elif s == 1:
                new_text.append(rules[item][1])
            else:
                new_text.append(rules[item][2])
        #Contexte sensitivity. If V is the last derivative produced 
        # probabilities are different
        elif item in rules and parsero[i+1]=='V':  # 0.4804	0.3513	0.1683
            s = np.random.uniform(low=0, high=1, size=None)
            if s < 0.4804:
                new_text.append(rules[item][0])
            elif 0.4804 <= s < 0.8317:
                new_text.append(rules[item][2])  
            else:    
                new_text.append(rules[item][1])     
        #Contexte sensitivity. If F is the last derivative produced 
        # probabilities are different
        elif item in rules and parsero[i+1]== 'F': # 0.8151	0.0015	0.1834
            s = np.random.uniform(low=0, high=1, size=None)
            if s < 0.8151:   
                new_text.append(rules[item][1])
            elif 0.8151 <= s < 0.8166:
                new_text.append(rules[item][0])  
            else:    
                new_text.append(rules[item][2])
        #Contexte sensitivity. If P is the last derivative produced 
        # probabilities are different        
        elif item in rules and parsero[i+1]== 'P': #0.0926	0.0418	0.8656
            s = np.random.uniform(low=0, high=1, size=None)
            if s < 0.0926 :
                new_text.append(rules[item][2])
            elif 0.0926 <= s < 0.1344 :
                new_text.append(rules[item][0])  
            else:    
                new_text.append(rules[item][1])  
                
        elif item in rules:
            new_text.append(rules[item])
        
        else:
            new_text.append(item)
    return ''.join(new_text)

def clsystemm(axioms, rules, iteration):
    x = axioms
    for i in range(iteration):
        x=iteratem(x, rules)
    return ''.join(x)

##Apply the previuos functions to get an array of files
#in which cells had equal probability of formation each cycle
lengths_array = np.random.randint(10,150,120)
random_files = []
for i in lengths_array:
    random_files.append(clsystemm('AB',rules_ctsens,i))
#Remove A and B
for x, y in enumerate(random_files): 
    random_files[x] = y.replace('B', '') 
for x, y in enumerate(random_files):     
    random_files[x] = y.replace('A', '') 
    
with open('../Data/Cell_files_data/contextmesicLsystem_edited_cells.txt', 'w') as output:
    for row in random_files:
        output.write(str(row) + '\n')

####Apply lsystem that generates Ray cells
rulesR = { "A" : ["AB","Z"] , "B" : ["V","F","P"], "Z": "ZR"}
##
def iteratex(axioms, rules):
    new_text = []
    parsero = list(axioms)
    parsero.append('')
    for i, item in enumerate(parsero):
        if item == '':
            pass
        #If string is just AB it will generate V,F or P with equal probability
        elif item in rules and parsero[i]=='A':
            s = random.uniform(0, 1)
            if s < 0.998:
                new_text.append(rules[item][0])
            else:
                new_text.append(rules[item][1])
        #   
        elif item in rules and parsero[i+1]=='':
            if parsero[i] == 'R':
                new_text.append(rules[item])
            elif parsero[i] == 'B':
                s = random.randint(0,3)
                if s == 0:
                    new_text.append(rules[item][0])
                elif s == 1:
                    new_text.append(rules[item][1])
                else:
                    new_text.append(rules[item][2])
            elif parsero[i] == 'Z':
                new_text.append(rules[item])
        elif item in rules and parsero[i]=='B':
            s = random.randint(0,3)
            if s == 0:
                new_text.append(rules[item][0])
            elif s == 1:
                new_text.append(rules[item][1])
            else:
                new_text.append(rules[item][2])
        elif item in rules and parsero[i]=='R':
            new_text.append(rules[item])
            
        elif item in rules:
            new_text.append(rules[item])
        
        else:
            new_text.append(item)
    return ''.join(new_text)

###
def lsystemr(axioms, rules, iteration):
    x = axioms
    for i in range(iteration):
        x=iteratex(x, rules)
    return ''.join(x)

#23% mean of files are ray cells
lengths_array = np.random.randint(10,150,120)
random_files = []
for i in lengths_array:
    s = random.uniform(0, 1)
    if s < 0.77:
        random_files.append(lsystemr('A',rulesR,i))
    else:
        random_files.append(lsystemr('Z',rulesR,i))

#Remove A,B and Z
for x, y in enumerate(random_files): 
    random_files[x] = y.replace('B', '') 
for x, y in enumerate(random_files):     
    random_files[x] = y.replace('A', '') 
for x, y in enumerate(random_files):     
    random_files[x] = y.replace('Z', '')
 
with open('../Data/Cell_files_data/Ray_Lsystem_edited_cells.txt', 'w') as output:
    for row in random_files:
        output.write(str(row) + '\n')


