#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 15:15:04 2018

@author: pedrolauand
"""
import cvxpy as cvx
import numpy as np

B=cvx.Variable((3,3),PSD=True) #Matriz do Elipsoide
d=cvx.Variable(3,1) #Centro do Elipsoide
constrains=[] #vinculos
n=int(input("Numero de linhas:"))#Numero de linhas da matriz do Politopo
c1,c2,c3=input("Primeira coluna:"),input("Segunda coluna:"),input("Terceira coluna:")#Colunas da matriz do Politopo
c1,c2,c3=c1.split(),c2.split(),c3.split()
for i in range(n):
    c1[i],c2[i],c3[i]=float(c1[i]),float(c2[i]),float(c3[i])
C=np.array([c1,c2,c3])
C=np.transpose(C)

b=input("Vetor independente:")#Vetor do Politopo(Ax-b<=0)
b=b.split()
for i in range(n):
    b[i]=float(b[i])

for i in range(n):
    C_t=np.array([[C[i][0]],[C[i][1]],[C[i][2]]])
    constrains.append(cvx.norm(cvx.matmul(B,C_t)) + cvx.matmul(C[i],d)<= b[i])
     #Vinculos do Problema
objective=cvx.Minimize(-1*(cvx.log_det(B))) #Funcao Objetiva

prob=cvx.Problem(objective,constrains)
prob.solve()
d=d.value
B=B.value
#vamos comecar o segundo passo
R1,D,R2=np.linalg.svd(B, full_matrices=True)
#Iremos encontrar as conjugações unitárias a partir de R1 e R2
#primeiro para R1
w1,v1=np.linalg.eig(R1)

for j in range(len(w1)):
    if w1[j].imag == 0 :
        v = v1[j]

trR1=0
for m in range(len(w1)):
    trR1=trR1 + w1[m]

alpha1 = np.arccos((trR1 - 1)*0.5)
#definindo as matrizes de pauli
sig1,sig2,sig3 = np.array([[0,1],[1,0]]),np.array([[0,-j],[j,0]]),np.array([[1,0],[0,-1]])
#como ja temos alpha e n podemos fazer a unitaria
U1=np.cos(alpha1*0.5)*np.identity(2, dtype=float) + j*np.sin(alpha1*0.5)*(v[0]*sig1 + v[1]*sig2 +v[2]*sig3)

w2,v2=np.linalg.eig(R2)


for j in range(len(w2)):
    if w2[j].imag == 0:
         v_p = v2[j]


trR2=0
for l in range(len(w2)):
    trR2=trR2 + w2[l]

alpha2 = np.arccos((trR2 - 1)*0.5)
U2=np.cos(alpha2*0.5)*np.identity(2, dtype=float) + j*np.sin(alpha2*0.5)*(v_p[0]*sig1 + v_p[1]*sig2 +v_p[2]*sig3)