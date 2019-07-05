#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 16:53:46 2019

@author: pedrolauand
"""

import numpy as np
###cenário(k,m,n) k partes m perguntas e n respostas para cada pergunta
####4plets (abxy) perguntas x e y respostas a e b
#(1111)-1:(2111)-2:(1121)-3:(2121)-4:(1211)-5:(2211)-6:(1221)-7:(2221)-8:(1112)-9
#(2112)-10:(1122)-11:(2122)-12:(1212)-13:(2212)-14:(1222)-15:(2222)-16
partes=int(input("Número de partes:"))
cenario=[]
for i in range(partes):
    cenario.append(input())
    cenario[i]=cenario[i].split()
    for j in range(len(cenario[i])):
        cenario[i][j]=int( cenario[i][j])
        
alpha=input("alfa :")
alpha=alpha.split()
for i in range(len(alpha)):
    alpha[i]=float(alpha[i])
#numero de perguntas len(cenario[j])
#numero de respostas sum(cenario[j])
listaB=[]
for m in range(len(cenario)):  
#construindo a base simetrica
    
    for i in range(sum(cenario[m])):
        if i==0:
            mhi=[[1/len(cenario[m])]]
        else:
            mhi[0].append(1/len(cenario[m]))
    nhi=[]
    for i in range(len(cenario[m])-1):
        nhi.append([0])
    
    for x in range(len(cenario[m])):
        for a in range(cenario[m][x]):
            for e in range(len(cenario[m])-1):
                if a==0 and x==0:
                    if x==e:
                        nhi[e]=[1]
                    elif x==e+1 :
                        nhi[e]=[-1]
                    else:
                        nhi[e]=[0]
                else:
                    if x==e:
                        nhi[e].append(1)
                    elif x==e+1 :
                        nhi[e].append(-1)
                    else:
                        nhi[e].append(0) 
    l=[]
    for e in range(len(cenario[m])):
        for i in range(cenario[m][e]-1):
            l.append([0])

    k=-1
    for e in range(len(cenario[m])):
        for c in range(cenario[m][e]-1):
            k=k+1     
            for x in range(len(cenario[m])):
                for a in range(cenario[m][x]): 
                    if a==0 and x==0:
                        if x!=e :
                            l[k]=[0]
                        elif a==c :
                            l[k]=[1]
                        elif a==c+1 :
                            l[k]=[-1]
                        else:
                            l[k]=[0]
                    else:
                        if x!=e :
                            l[k].append(0)
                        elif a==c :
                            l[k].append(1)
                        elif a==c+1 :
                            l[k].append(-1)
                        else:
                            l[k].append(0)
    base=mhi+l+nhi
    listaB.append(base)

gamma,G=[],[]
for i in listaB[1]:
    for j in listaB[0]:
        gamma.append(np.kron(i,j))
for i in range(len(alpha)):
    G.append(np.dot(gamma[i],alpha)/(np.linalg.norm(gamma[i]))**2)
print(gamma[3],gamma[7],gamma[11],gamma[12],gamma[13],gamma[14],gamma[15])

             
    

