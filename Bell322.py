
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 16:53:46 2019

@author: pedrolauand
"""
## Esse código acha o mínimo lexicográfico já sem as degenerescências criadas pelo subspaço de normalização e não sinalizante para uma dada
## desigualdade do cenário de Bell (3,2,2).(tripartido homogêneo)
## A entrada deve ser feita somente pelo vetor de coeficientes que acompanham as probabilidades da desigualdade
## Essa entrada deve ser posta como texto sem indicar o sinal de mais, mas indicando o sinal de menos na frente do coeficiente sem espaco entre o "-" e o número
##
import numpy as np
from sympy.combinatorics import Permutation
from sympy.combinatorics.perm_groups import PermutationGroup
import copy

partes=3
cenario=[[2,2],[2,2]]

alpha=input("alfa :")#vetor de coeficientes
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
print(Gamma,"<=",beta)##Desigualdade após a projeção no subespaço não sinalizante

##Construindo o grupo de todas as permutações do cenário

p1=Permutation([0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15,16,20,24,28,17,21,25,29,18,22,26,30,19,23,27,31,32,36,40,44,33,37,41,45,34,38,42,46,35,39,43,47,48,52,56,60,49,53,57,61,50,54,58,62,51,55,59,63])
p2=Permutation([2,3,0,1,6,7,4,5,10,11,8,9,14,15,12,13,18,19,16,17,22,23,20,21,26,27,24,25,30,31,28,29,34,35,32,33,38,39,36,37,42,43,40,41,46,47,44,45,50,51,48,49,54,55,52,53,58,59,56,57,62,63,60,61])
p3=Permutation([1,0,2,3,5,4,6,7,9,8,10,11,13,12,14,15,17,16,18,19,20,21,22,23,25,24,26,27,29,28,30,31,33,32,34,35,37,36,38,39,41,40,42,43,45,44,46,47,49,48,50,51,53,52,54,55,57,56,58,59,61,60,62,63])
p4=Permutation([0,1,3,2,4,5,7,6,8,9,11,10,12,13,15,14,16,17,19,18,20,21,23,22,24,25,27,26,28,29,31,30,32,33,35,34,36,37,39,38,40,41,43,42,44,45,47,46,48,49,51,50,52,53,55,54,56,57,59,58,60,61,63,62])
p5=Permutation([8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7,24,25,26,27,28,29,30,31,16,17,18,19,20,21,22,23,40,41,42,43,44,45,46,47,32,33,34,35,36,37,38,39,56,57,58,59,60,61,62,63,48,49,50,51,52,53,54,55])
p6=Permutation([4,5,6,7,0,1,2,3,8,9,10,11,12,13,14,15,20,21,22,23,16,17,18,19,24,25,26,27,28,29,30,31,36,37,38,39,32,33,34,35,40,41,42,43,44,45,46,47,52,53,54,55,48,49,50,51,56,57,58,59,60,61,62,63])
p7=Permutation([0,1,2,3,4,5,6,7,12,13,14,15,8,9,10,11,16,17,18,19,20,21,22,23,28,29,30,31,24,25,26,27,32,33,34,35,36,37,38,39,44,45,46,47,40,41,42,43,48,49,50,51,52,53,54,55,60,61,62,63,56,57,58,59])
p8=Permutation([0,1,2,3,16,17,18,19,32,33,34,35,48,49,50,51,4,5,6,7,20,21,22,23,36,37,38,39,52,53,54,55,8,9,10,11,24,25,26,27,40,42,43,43,56,57,58,59,12,13,14,15,28,29,30,31,44,45,46,47,60,61,62,63])
p9=Permutation([32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31])
p10=Permutation([16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63])
p11=Permutation([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47])
G=PermutationGroup([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11])
print(G.order())
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
for j in range(len(Gamma)-1): ##Achando o mínimo lexicográfico
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
print(alpha)#Representante da órbita
