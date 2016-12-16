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
import DFR_1D

### Parameters

# SD order
p       =  3

#Solver
# method=0 for Convection  !!!! c=/0!!!
# method=1 for Diffusion
#method=2 for Convection+Diffusion

#method=    0 -> Pas besoin de ça (AS)


# Stability criterion : CFL 
CFL=       0.7

#Final time
Tfin   = 0.002 


# Velocity c (m/s) and diffusion D (m^2/s)
c      =  10.
D      = 0.

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
bcond = 1;
# If Dirichlet conditions specify values
yL=0
yR=0

### Computating solution
x0, sol0, x, sol, niter = DFR_1D.main(p,CFL,Tfin,c,D,init,grad_init,bcond,yL*(1-bcond),yR*(1-bcond))



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
    if init==0:
        yth[i]=1/m.sqrt((4*D*Tfin/0.0001)+1)*m.exp(-(x[i]-c*Tfin+0.1)**2/(4*D*Tfin+0.0001))
    elif init==1:
        if method==0:
            yth[i]=(1-m.erf((x[i]-c*(Tfin+grad_init))))/2
        else:
            yth[i]=(1-m.erf((x[i]-c*(Tfin+grad_init))/(2*m.sqrt(D*(Tfin+grad_init)))))/2
    
    file.write(str(x[i])+",")
    file.write(str(yth[i])+",")
    file.write(str(sol[i]) + "\n")

file.close()


### Solution plotting
plt.plot(x0,sol0,'k-')
plt.plot(x,sol, 'r-')
plt.plot(x,yth, 'b-')
plt.legend(['Initial', 'Convected','Theoretical'])
plt.show()

