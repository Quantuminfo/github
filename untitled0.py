#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 16:53:46 2019

@author: pedrolauand
"""

###cenário(k,m,n) k partes m perguntas e n respostas para cada pergunta
####4plets (abxy) perguntas x e y respostas a e b
#(1111)-1:(2111)-2:(1121)-3:(2121)-4:(1211)-5:(2211)-6:(1221)-7:(2221)-8:(1112)-9
#(2112)-10:(1122)-11:(2122)-12:(1212)-13:(2212)-14:(1222)-15:(2222)-16
exp=[int(input()),int(input()),int(input())]
alpha=input()
alpha=alpha.split()
for i in range(len(alpha)):
    alpha[i]=float(alpha[i])
    
#construindo a base simetrica
mhi=[]
for i in range(exp[1]*exp[2]):
    mhi.append(1/exp[0])
nhi=[]
for x in range(exp[1]):
    for a in range(exp[2]):
        for e in range(exp[1]-1):
            if x==e:
                nhi.append(1)
            elif x==e+1 :
                nhi.append(-1)
print(nhi,mhi)            
    

