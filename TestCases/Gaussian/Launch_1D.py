#!/usr/bin/env python
# -*- coding: utf-8 -*-
''' author: R. Pile & A. Schmalzried '''
''' 19/11/2016 '''
''' Launch_1D.py : Script for Gaussian test case (with Gauss-Lobatto SPs) '''
''' This work follows python scripts developped by J. Vanharen and V. Joncquieres, under the supervision of G. Puigt '''


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

# Number of cycles
NCycles = 1000 

# Spectral method order
p = 4

# Number of cells
N = 5

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
plt.ylim([-0.5,2.5])
plt.plot(x,yth,'k')
plt.plot(x0SD,sol0SD,'d',color=[1,0,0],markersize = 10,markerfacecolor='None',markeredgecolor=[1,0,0],markeredgewidth=2)
plt.plot(x0FR,sol0FR,'+',color=[0.5,0,0.5],markersize = 10,markeredgewidth=2)
plt.plot(x0DFR,sol0DFR,'o',color=[0,0,1],markersize = 10,markerfacecolor='None',markeredgecolor=[0,0,1],markeredgewidth=2)
if(fou):
    plt.plot(x0FD,sol0FD,'x',color=[0,0.4,0.4],markersize = 10,markeredgewidth=2)
    plt.legend(['True solution','SD','FR','DFR','FOU'],loc='upper right')
else:
    plt.legend(['True solution','SD','FR','DFR'],loc='upper right')


plt.plot(xSD,solSD,color=[1,0,0])
plt.plot(xFR,solFR,color=[0.5,0,0.5])
plt.plot(xDFR,solDFR,color=[0,0,1])
if(fou):
    plt.plot(xFD,solFD,color=[0,0.4,0.4])




plt.figure()

Nticks = 20


NSD = int(np.floor((niterSD+1)/Nticks))
NFR = int(np.floor((niterFR+1)/Nticks))
NDFR = int(np.floor((niterDFR+1)/Nticks))
if fou:
    NFD = int(np.floor((niterFD+1)/Nticks))

VSD = np.linspace(0,niterSD+1,niterSD+2)
VFR = np.linspace(0,niterFR+1,niterFR+2)
VDFR = np.linspace(0,niterDFR+1,niterDFR+2)
if fou :
    VFD = np.linspace(0,niterFD+1,niterFD+2)



plt.plot(VSD[0:-2:NSD]/10000.,l2SD[0:-2:NSD],'d',color=[1,0,0],markersize = 10,markerfacecolor='None',markeredgecolor=[1,0,0],markeredgewidth=2)
plt.plot(VFR[0:-2:NFR]/10000.,l2FR[0:-2:NFR],'+',color=[0.5,0,0.5],markersize = 10,markeredgewidth=2)
plt.plot(VDFR[0:-2:NDFR]/10000.,l2DFR[0:-2:NDFR],'o',color=[0,0,1],markersize = 10,markerfacecolor='None',markeredgecolor=[0,0,1],markeredgewidth=2)
if(fou):
    plt.plot(VFD[0:NFD*Nticks:NFD],l2FD[0:-1:NFD],'x',color=[0,0.5,0.75],markersize = 10,markeredgewidth=2)
    plt.legend(['SD','FR','DFR','FOU'],loc='upper right')
else:
    plt.legend(['SD','FR','DFR'],loc='upper right')


plt.plot(np.linspace(0,niterSD,niterSD+1)/10000.,l2SD[0:-1],color=[1,0,0])
plt.plot(np.linspace(0,niterFR,niterFR+1)/10000.,l2FR[0:-1],color=[0.5,0,0.5])
plt.plot(np.linspace(0,niterDFR,niterDFR+1)/10000.,l2DFR[0:-1],color=[0,0,1])
if(fou):
    plt.plot(np.linspace(0,niterFD+1,niterFD+2),l2FD,color=[0,0.5,0.75])
plt.title("L2-norm")
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


