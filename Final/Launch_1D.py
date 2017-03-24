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
import FR_1D
import SDcombu_diff
import solutionPoints

### Parameters

# Scheme
# scheme = 'SD', 'FR', 'DFR'
scheme = 'SD'

# DFR order
p = 3

# Number of cells
N = 50  

# Type of SP
# SP = 1    --> Gauss-Lobatto
# SP = 2    --> Gauss
# SP = 3    --> SPs with z1 from Romero et al.
SPchoice = 2
z1 = 0.                 # z1 from Romero et al.

# Stability criterion : CFL 
CFL = 0.8

#Final time
Tfin = 0.05

# Domain Length
L = 1.0

# Velocity c (m/s) and diffusion D (m^2/s)
c = 10.
D = 0.

#Initialization 
# init='Gauss' --> Gaussian law 
#init='Constant' --> Constant law 
# init = 'Triangle' --> Linear Tipi law
#init='RectErf'--> Rectangular step with erf function
init='Sine'


# Gradient for initialize the erf 
grad_init=10**(-12)

#boundary conditions on the left side bcondL and the right side bcondR
#bcond*=0 for a Dirichlet BC
#bcond*=1 for a periodic BC
bcond = 1 
# If Dirichlet conditions specify values
yL=0.1
yR=0.


# Time integration scheme
# timeIntegration='RK4'     --> Fourth order Runge-Kutta scheme
# timeIntegration='RK6low'  --> Low storage sixth order Runge-Kutta scheme 
timeIntegration='RK6low'

#Cell spacing :
# cellmask = 'Regular' -> evenly spaced cells
# cellmask = 'Irregular' -> unevenly customisable cell spacing
cellmask = 'Regular'



    ########################################################
    #                                                      #
    #                        SD                            #
    #                                                      #
    ########################################################

#Solver (ONLY FOR SD)
# method=0 for Convection  !!!! c=/0!!!
# method=1 for Diffusion
#method=2 for Convection+Diffusion
method = 0

#boundary conditions on the left side bcondL and the right side bcondR
#bcond*=0 for a Dirichlet BC
#bcond*=1 for a periodic BC
bcondL=1
bcondR=1




    ########################################################
    #                                                      #
    #                        FR                            #
    #                                                      #
    ########################################################

# Correction function (ONLY FOR FR)
 # corFun = 0 : g2 correction function
 # corFun = 1 : gr and gl correction functions
 # ...
corFun = 0.


    ########################################################
    #                                                      #
    #                        DFR                           #
    #                                                      #
    ########################################################

 # Penalizing parameter (ONLY FOR DFR)
tau = 0.


#############################################################################################################
####################################    Computing solution      #############################################
#############################################################################################################
solutionPoints.writeSP(p,SPchoice,z1)
if scheme == 'DFR':
    x0, sol0, x, sol, niter, l2 = DFR_1D.main(p,CFL,Tfin,c,D,init,grad_init,bcond,yL,yR,N,L,tau,timeIntegration,cellmask)
elif scheme =='SD':
    x0, sol0, x, sol, niter, l2 = SDcombu_diff.main(p,method,CFL,Tfin,c,D,init,grad_init,bcondL,bcondR,yL,yR,N,L,timeIntegration,cellmask)
elif scheme == 'FR':
    x0, sol0, x, sol, niter, l2 = FR_1D.main(p,CFL,Tfin,c,D,init,grad_init,bcond,yL,yR,N,L,corFun,timeIntegration,cellmask)
else:
    print 'Unrecognized method'
    quit()

# Write initialization
file=open("results_init.txt","w")
for i in range(len(x0)):
    file.write(str(x0[i])+",")
    file.write(str(sol0[i]) + "\n")

file.close()

# Write results, need to modify filename 
file = open(scheme+"_results_CFL="+str(CFL)+"_p="+str(p)+".txt", "w")   
yth=np.zeros(len(x))
for i in range(len(x)):
    if init=='Gauss':
        if bcond == 0:
            yth[i]=1/m.sqrt((4*D*20*Tfin)+1)*m.exp(-(x[i]-c*Tfin)**2/(4*D*Tfin+0.05))
        else:
            yth[i]=1/m.sqrt((4*D*20*Tfin)+1)*m.exp(-(np.mod(x[i]-c*Tfin+L/2,L)-L/2)**2/(4*D*Tfin+0.05))
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


