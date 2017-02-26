#!/usr/bin/env python
# -*- coding: utf-8 -*-
''' author: R. Pile & A. Schmalzried '''
''' 19/11/2016 '''
''' Launch_1D.py : Script for DFR method '''
import sys
import numpy as np
import matplotlib.pyplot as plt
import math as m
from scipy import integrate
import SDcombu_diff
import FR_1D
import DFR_1D
import FD_1D
import time

### Parameters

# Use First Order Upwind ? 
fou = False

# Plot from period nstart to period nend
nstart = 3.
nend = 5. 

# Sin wave number
ksin = np.pi/4

# Number of cycles
NCycles = 5.

# Spectral method order
p = 4

# Number of cells
N = 10

# Stability criterion : CFL 
CFL = 0.8                                  # La condition CFL semble plus restrictive que CFL < 1 (GR)

# Domain Length
L = 20.0

# Velocity c (m/s) and diffusion D (m^2/s)
c = 1.
D = 0. #2.5e-3    


#Final time
Tfin = L/c*NCycles

#Solver
# method=0 for Convection  !!!! c=/0!!!
# method=1 for Diffusion
#method=2 for Convection+Diffusion

method=    0

 # Penalizing parameter
tau = 0.

 # Correction function
 # corFun = 0 : g2 correction function
 # corFun = 1 : gr and gl correction functions
 # ...
corFun = 1.


#Initialization 
# init='Gauss' --> Gaussian law 
#init='Constant' --> Constant law 
# init = 'Triangle' --> Linear Tipi law
#init='RectErf'--> Rectangular step with erf function

init='Constant'


# Gradient for initialize the erf 
grad_init=10**(-12)

#boundary conditions on the left side bcondL and the right side bcondR
#bcond*=0 for a Dirichlet BC
#bcond*=1 for a periodic BC
bcond = 0.
# If Dirichlet conditions specify values
yL=0.
yR=0.


# Time integration scheme
# timeIntegration='RK4'     --> Fourth order Runge-Kutta scheme
# timeIntegration='RK6low'  --> Low storage sixth order Runge-Kutta scheme 
timeIntegration='RK4'

#Cell spacing :
# cellmask = 'Regular' -> evenly spaced cells
# cellmask = 'Irregular' -> unevenly customisable cell spacing
cellmask = 'Regular'

### Computing solution
t0 = time.time()
x0SD, sol0SD, xSD, solSD, niterSD, l2SD = SDcombu_diff.main(p,method,CFL,Tfin,c,D,0,grad_init,0.0,1.0,yL,yR,N,L,timeIntegration,cellmask,ksin)
tSD = time.time()-t0


t0 = time.time()
x0FR, sol0FR, xFR, solFR, niterFR, l2FR = FR_1D.main(p,CFL,Tfin,c,D,init,grad_init,bcond,yL,yR,N,L,corFun,timeIntegration,cellmask,ksin)
tFR = time.time()-t0

t0 = time.time()
x0DFR, sol0DFR, xDFR, solDFR, niterDFR, l2DFR = DFR_1D.main(p,CFL,Tfin,c,D,init,grad_init,bcond,yL,yR,N,L,tau,timeIntegration,cellmask,ksin)
tDFR = time.time()-t0

if(fou):
    t0 = time.time()
    x0FD, sol0FD, xFD, solFD, niterFD, l2FD = FD_1D.main(p,CFL,Tfin,c,D,init,grad_init,bcond,yL,yR,N,L,timeIntegration,cellmask,ksin)
    tFD = time.time()-t0



x = np.linspace(0.,L,10*N*p)

# Write results, need to modify filename 
file = open("results_CFL="+str(CFL)+"_p="+str(p)+".txt", "w")   
yth=np.zeros(len(x))
for i in range(len(x)):
    if init=='Gauss':
        if bcond == 0:
            yth[i]=1/m.sqrt((4*D*Tfin/0.0001)+1)*m.exp(-(x[i]-c*Tfin+0.1)**2/(4*D*Tfin+0.0001))
        else:
            yth[i]=m.exp(-20*(np.mod(x[i]-c*Tfin+L/2,L)-L/2)**2)
    elif init=='Constant':
        yth[i] = 10. 
    elif init=='RectErft':
        if method==0:
            yth[i]=(1-m.erf((x[i]-c*(Tfin+grad_init))))/2
        else:
            yth[i]=(1-m.erf((x[i]-c*(Tfin+grad_init))/(2*m.sqrt(D*(Tfin+grad_init)))))/2
    elif init=='Sine':
            yth[i]=m.cos(2*2*np.pi*(x[i]-c*Tfin)) 


file.close()

#print 'Solution initiale'
#print sol0
#print 'Solution calculee'
#print sol

#for i in range(len(x)):
#    print sol0[i]

### Solution plotting

plt.rcParams.update({'font.size': 18})

plt.figure()
plt.ylim([-2,2])
plt.plot(xSD,solSD)
plt.plot(xFR,solFR)
plt.plot(xDFR,solDFR)
if(fou):
    plt.plot(xFD,solFD)
    plt.legend(['SD','FR','DFR','FOU'],loc='upper right')
else:
    plt.legend(['SD','FR','DFR'],loc='upper right')
plt.title("Solution at ct/L = "+str(NCycles))
plt.xlabel("x")
plt.ylabel("u")

plt.figure()
plt.plot(np.linspace(0,niterSD+1,niterSD+2),l2SD)
plt.plot(np.linspace(0,niterFR+1,niterFR+2),l2FR)
plt.plot(np.linspace(0,niterDFR+1,niterDFR+2),l2DFR)
if(fou):
    plt.plot(np.linspace(0,niterFD+1,niterFD+2),l2FD)
    plt.legend(['SD','FR','DFR','FOU'],loc='upper right')
else:
    plt.legend(['SD','FR','DFR'],loc='upper right')
plt.title("Integral difference")
plt.xlabel("Iteration")
plt.ylabel("L2-norm")
xlim1 = int(nstart*(niterSD+2)/float(NCycles))
xlim2 = int(nend*(niterSD+2)/float(NCycles))
ylim1 = min([min(l2SD[xlim1:xlim2]),min(l2FR[xlim1:xlim2]),min(l2DFR[xlim1:xlim2])])
ylim2 = max([max(l2SD[xlim1:xlim2]),max(l2FR[xlim1:xlim2]),max(l2DFR[xlim1:xlim2])])
plt.xlim([xlim1,xlim2])
plt.ylim([1.5*ylim1,1.5*ylim2])

# Write initialization
file=open("times.txt","w")
file.write(str(tSD) + "\n")
file.write(str(tFR) + "\n")
file.write(str(tDFR) + "\n")
if(fou):
    file.write(str(tFD) + "\n")
file.close()


plt.show()


