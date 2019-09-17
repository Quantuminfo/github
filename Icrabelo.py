"""
Created on Fri Oct 12 15:15:04 2018

@author: pedrolauand
"""
import cvxpy as cvx
import numpy as np


B=cvx.Variable((3,3),PSD=True) #Matriz do Elipsoide
d=cvx.Variable(3,1) #Centro do Elipsoide
constrains=[] #vinculos

inputiq= open("inputiq.txt")
q=inputiq.readlines()
inputiq.close()

c1,c2,c3=q[0],q[1],q[2]#Colunas da matriz do Politopo
c1,c2,c3=c1.split(),c2.split(),c3.split()
C=np.array([c1,c2,c3])
C=np.transpose(C)

for i in range(len(C)):
    c1[i],c2[i],c3[i]=float(c1[i]),float(c2[i]),float(c3[i])



for i in range(len(C)):
    c1[i],c2[i],c3[i]=float(c1[i]),float(c2[i]),float(c3[i])


b=q[3] #Vetor do Politopo(Ax-b<=0)
b=b.split()
for i in range(len(C)):
    b[i]=float(b[i])

for i in range(len(C)):
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

trR1=R1[1][1]+R1[2][2]+R1[0][0]

alpha1 = np.arccos((trR1 - 1)*0.5)
#definindo as matrizes de pauli
sig1,sig2,sig3 = np.array([[0,1],[1,0]]),np.array([[0,-j],[j,0]]),np.array([[1,0],[0,-1]])
#como ja temos alpha e n podemos fazer a unitaria
U1=np.cos(alpha1*0.5)*np.identity(2, dtype=float) + j*np.sin(alpha1*0.5)*(v[0]*sig1 + v[1]*sig2 +v[2]*sig3)

w2,v2=np.linalg.eig(R2)


for j in range(len(w2)):
    if w2[j].imag == 0:
         v_p = v2[j]


trR2=R2[1][1]+R2[2][2]+R2[0][0]


alpha2 = np.arccos((trR2 - 1)*0.5)
U2=np.cos(alpha2*0.5)*np.identity(2, dtype=float) + j*np.sin(alpha2*0.5)*(v_p[0]*sig1 + v_p[1]*sig2 +v_p[2]*sig3)
#agora so precisamos dos operadores de Kraus
choi = 0.5*(np.array([[(1+D[2]+d[2]),(d[0]+j*d[1]),0,(D[0]+D[1])],[(d[0]-j*d[1]),(1-D[2]-d[2]),(D[0]-D[1]),0],[0,(D[0]-D[1]),(1-D[2]+d[2]),(d[0]+j*d[1])],[(D[0]+D[1]),0,(d[0]-j*d[1]),(1+D[2]-d[2])]]))
w3,v3=np.linalg.eig(choi)
#definindo os projetores
P1,P2 = np.array([[1,0,0,0],[0,1,0,0]]),np.array([[0,0,1,0],[0,0,0,1]])
P1,P2=np.reshape(P1,(2,4)),np.reshape(P2,(2,4))
v3=np.transpose(v3)
#definindo K(op à esquerda)
K=[0,1,2,3]


for j in range(len(w3)):
    if w3[j] <= 10**(-3):
        w3[j]=0
    for i in range(len(w3)):
        if abs(v3[j][i]) <= 10**(-3):
            v3[j][i]=0
    v3[j]=np.transpose(v3[j])
    K[j]=np.array([[np.matmul(P1,v3[j])],[np.matmul(P2,v3[j])]])
    K[j]=np.transpose(K[j])
    K[j]=np.sqrt(w3[j])*K[j]
    K[j]=np.reshape(K[j],(2,2))

#operadores de Kraus à direita K_dagger
K_dagger=[0,1,2,3]

for j in range(len(K)):
    K_dagger[j]=np.matrix.conjugate(K[j])
    K_dagger[j]=np.transpose(K_dagger[j])
    K_dagger[j]=np.reshape(K_dagger[j],(2,2))

ro=np.array([[1,0],[0,0]])
ro=np.reshape(ro,(2,2))
au=[0,0,0,0]
aux=[0,0,0,0]
for i in range(len(K)):
    au[i]=np.matmul(K[i],ro)
    aux[i]=np.matmul(au[i],K_dagger[i])
rop=aux[0]+aux[1]+aux[2]+aux[3]
I=[0,0,0,0]
for i in range(len(K)):
    I[i]=np.matmul(K_dagger[i],K[i])
I=I[0]+I[1]+I[2]+I[3]
print(rop)
print(I)
