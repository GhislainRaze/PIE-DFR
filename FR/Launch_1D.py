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
import FR_1D

### Parameters

# FR order
p = 5

#Number of cells
N = 100

# Stability criterion : CFL 
CFL = 0.8

#Final time
Tfin = 0.1

# Domain Length
L = 1.0

# Velocity c (m/s) and diffusion D (m^2/s)
c = 10.
D = 2.5e-3

 # Correction function
 # corFun = 0 : g2 correction function
 # corFun = 1 : gr and gl correction functions
 # ...
corFun = 1.

#Initialization 
# init = 'Gauss'    --> Gaussian law 
# init = 'Constant' --> Constant law 
# init = 'Triangle' --> Linear Tipi law
# init = 'RectErf'  --> Rectangular step with erf function
# init = 'Sine'     --> Sine
init='Gauss'


# Gradient for initialize the erf 
grad_init=10**(-12)

#boundary conditions on the left side bcondL and the right side bcondR
#bcond*=0 for a Dirichlet BC
#bcond*=1 for a periodic BC
bcond = 1


# If Dirichlet conditions specify values
yL=0.5
yR=0.

# Time integration scheme
# timeIntegration='RK4'     --> Fourth order Runge-Kutta scheme
# timeIntegration='RK6low'  --> Low storage sixth order Runge-Kutta scheme 
timeIntegration='RK6low'

#Cell spacing :
# cellmask = 'Regular'  --> evenly spaced cells
# cellmask = 'Irregular'--> unevenly customisable cell spacing
cellmask = 'Regular'

### Computing solution
x0, sol0, x, sol, niter = FR_1D.main(p,CFL,Tfin,c,D,init,grad_init,bcond,yL,yR,N,L,corFun,timeIntegration,cellmask)



# Write initialization
file=open("results_init.txt","w")
for i in range(len(x0)):
    file.write(str(x0[i])+",")
    file.write(str(sol0[i]) + "\n")

file.close()

# Write results, need to modify filename 
file = open("results_CFL="+str(CFL)+"_p="+str(p)+".txt", "w")   
yth=np.zeros(len(x))
for i in range(len(x)):
    if init=='Gauss':
        if bcond == 0:
            yth[i]=1/m.sqrt((4*D*Tfin/0.0001)+1)*m.exp(-(x[i]-c*Tfin+0.1)**2/(4*D*Tfin+0.0001))
        else:
            yth[i]=1/m.sqrt((4*D*Tfin/0.0001)+1)*m.exp(-(np.mod(x[i]-c*Tfin+0.1+L/2,L)-L/2)**2/(4*D*Tfin+0.0001))
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

### Solution plotting
plt.plot(x0,sol0,'k-')
plt.plot(x,sol, 'r-')
plt.plot(x,yth, 'b--')
plt.legend(['Initial', 'Convected','Theoretical'])
plt.show()


