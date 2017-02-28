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
import solutionPoints

### Parameters

# Use First Order Upwind ? 
fou = False

# Number of cycles
NCycles = 10

# Spectral method order
p = 4

# Number of cells
N = 5

# Type of SP
# SP = 1    --> Gauss-Lobatto
# SP = 2    --> Gauss
# SP = 3
SPchoice = 1

solutionPoints.writeSP(p,SPchoice)


# Stability criterion : CFL 
CFL = 0.8                                  # La condition CFL semble plus restrictive que CFL < 1 (GR)

# Domain Length
L = 2.0

# Velocity c (m/s) and diffusion D (m^2/s)
c = 10.
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

init='Gauss'


# Gradient for initialize the erf 
grad_init=10**(-12)

#boundary conditions on the left side bcondL and the right side bcondR
#bcond*=0 for a Dirichlet BC
#bcond*=1 for a periodic BC
bcond = 1
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
x0SD, sol0SD, xSD, solSD, niterSD, l2SD = SDcombu_diff.main(p,method,CFL,Tfin,c,D,0,grad_init,bcond,bcond,yL,yR,N,L,timeIntegration,cellmask)
tSD = time.time()-t0


t0 = time.time()
x0FR, sol0FR, xFR, solFR, niterFR, l2FR = FR_1D.main(p,CFL,Tfin,c,D,init,grad_init,bcond,yL,yR,N,L,corFun,timeIntegration,cellmask)
tFR = time.time()-t0

t0 = time.time()
x0DFR, sol0DFR, xDFR, solDFR, niterDFR, l2DFR = DFR_1D.main(p,CFL,Tfin,c,D,init,grad_init,bcond,yL,yR,N,L,tau,timeIntegration,cellmask)
tDFR = time.time()-t0

if(fou):
    t0 = time.time()
    x0FD, sol0FD, xFD, solFD, niterFD, l2FD = FD_1D.main(p,CFL,Tfin,c,D,init,grad_init,bcond,yL,yR,N,L,timeIntegration,cellmask)
    tFD = time.time()-t0



x = np.linspace(-0.5*L,0.5*L,10*N*p)

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
plt.plot(x,yth, 'k-')
plt.plot(xSD,solSD)
plt.plot(xFR,solFR)
plt.plot(xDFR,solDFR)
if(fou):
    plt.plot(xFD,solFD)
    plt.legend(['Initial','SD','FR','DFR','FOU'],loc=0)
else:
    plt.legend(['Initial','SD','FR','DFR'],loc=0)
plt.title("Solution")
plt.xlabel("x")
plt.ylabel("u")

plt.figure()
plt.plot(np.linspace(0,niterSD+1,niterSD+2)/10000.,l2SD)
plt.plot(np.linspace(0,niterFR+1,niterFR+2)/10000.,l2FR)
plt.plot(np.linspace(0,niterDFR+1,niterDFR+2)/10000.,l2DFR)
if(fou):
    plt.plot(np.linspace(0,niterFD+1,niterFD+2),l2FD)
    plt.legend(['SD','FR','DFR','FOU'],loc=0)
else:
    plt.legend(['SD','FR','DFR'],loc=0)
plt.title("L2-norm of the solution")
plt.xlabel("Iteration (10^4)")
plt.ylabel("L2-norm")

# Write initialization
file=open("times.txt","w")
file.write(str(tSD) + "\n")
file.write(str(tFR) + "\n")
file.write(str(tDFR) + "\n")
if(fou):
    file.write(str(tFD) + "\n")
file.close()


plt.show()


