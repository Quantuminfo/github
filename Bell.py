
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 16:53:46 2019

@author: pedrolauand
"""

import numpy as np
from sympy.combinatorics import Permutation
from sympy.combinatorics.perm_groups import PermutationGroup
import copy


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
projlistaB=[]
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
    aux=copy.deepcopy(nhi)

    for i in range(len(aux)):
        for j in range(len(aux[i])):
            aux[i][j]=0
    projbase=mhi+l+aux
    projlistaB.append(projbase)
projauxiliar=[projlistaB[partes-1]]
for m in range(partes-1):
    projauxiliar1=[]
    for i in projauxiliar[m]:
        for j in projlistaB[partes-2-m]:
            projauxiliar1.append(np.kron(i,j))
    projauxiliar.append(projauxiliar1)
projgamma=projauxiliar[partes-1] 

auxiliar=[listaB[partes-1]]
for m in range(partes-1):
    auxiliar1=[]
    for i in auxiliar[m]:
        for j in listaB[partes-2-m]:
            auxiliar1.append(np.kron(i,j))
    auxiliar.append(auxiliar1)
gamma=auxiliar[partes-1]
gamma=np.transpose(gamma)
gamma=np.linalg.inv(gamma)
Gamma=np.matmul(gamma,alpha)
for i in range(len(gamma)):
    if all(elem==0 for elem in projgamma[i]):
        Gamma[i]=0
beta=-Gamma[0]
Gamma[0]=0
print(Gamma,"<=",beta)

list=[]
for k in range(len(cenario)):
    list.append(k)
for k in range(len(cenario)):
    list[k]=[]
    for x in range(len(cenario[k])):
        for a in range(cenario[k][x]):
            list[k].append([a,x])
print(list)
  
p1=Permutation(1, 4)(2, 8)(3, 12)(6, 9)(7, 13)(11, 14)
p2=Permutation(0, 2)(1, 3)(4, 6)(5, 7)(8, 10)(9, 11)(12, 14)(13, 15)
p3=Permutation(0, 8)(1, 9)(2, 10)(3, 11)(4, 12)(5, 13)(6, 14)(7, 15)
p4=Permutation(0, 1)(4, 5)(8, 9)(12, 13)
p5=Permutation(2, 3)(6, 7)(10, 11)(14, 15)
p6=Permutation(0, 4)(1, 5)(2, 6)(3, 7)
p7=Permutation(8, 12)(9, 13)(10, 14)(11, 15)
G=PermutationGroup(p1,p2,p3,p4,p5,p6,p7)
lista=[]
Stab=[G]
for j in range(len(alpha)):
    lista.append(j)
    Stab.append(G.pointwise_stabilizer(lista))
    
U=[]
for i in range(len(Stab)-1):
    U.append(Stab[i].coset_transversal(Stab[i+1]))
H=[[G.identity]]
for i in range(len(Gamma)-1):
    H.append([])
for j in range(len(Gamma)-1):
    m=np.infty
    for u in U[j]:
        for h in H[j]:
            perm=[i^u^h for i in range(len(Gamma))]
            if Gamma[perm[j]]<m:
                m=Gamma[perm[j]]
                H[j+1]=[]
            if Gamma[perm[j]]==m:
                H[j+1].append(u*h)
permu=[i^H[len(Gamma)-1][0] for i in range(len(Gamma))]
aux1=copy.deepcopy(Gamma)
for i in range(len(Gamma)):
    alpha[i]=aux1[permu[i]]
print(alpha)
