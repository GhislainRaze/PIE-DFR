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
NCycles = 1000

# Spectral method order
p = 5

# Number of cells
N = 5

# Type of SP
# SP = 1    --> Gauss-Lobatto
# SP = 2    --> Gauss
# SP = 3
SPchoice = 3
if p == 3:
    z1SP = 0.339842589774454
    z1DG = np.sqrt(3./7.-2./7.*np.sqrt(6./5.))
    z11 = 0.3402
    z12 = 0.35
elif p == 4:
    z1SP = 0.538323058771738
    z1DG = 1./3.*np.sqrt(5.-2.*np.sqrt(10./7.))
    z11 = 0.5387
    z12 = 0.55
elif p == 5:
    z1SP = 0.238431046729096
    z1DG = 0.238619186083197
    z11 = 0.2389
    z12 = 0.25



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
solutionPoints.writeSP(p,SPchoice,z1SP)
t0 = time.time()
x0DFR1, sol0DFR1, xDFR1, solDFR1, niterDFR1, l2DFR1 = DFR_1D.main(p,CFL,Tfin,c,D,init,grad_init,bcond,yL,yR,N,L,tau,timeIntegration,cellmask)
tDFR1 = time.time()-t0

solutionPoints.writeSP(p,SPchoice,z1DG)
t0 = time.time()
x0DFR1, sol0DFR2, xDFR2, solDFR2, niterDFR2, l2DFR2 = DFR_1D.main(p,CFL,Tfin,c,D,init,grad_init,bcond,yL,yR,N,L,tau,timeIntegration,cellmask)
tDFR2 = time.time()-t0

solutionPoints.writeSP(p,SPchoice,z11)
t0 = time.time()
x0DFR1, sol0DFR3, xDFR3, solDFR3, niterDFR3, l2DFR3 = DFR_1D.main(p,CFL,Tfin,c,D,init,grad_init,bcond,yL,yR,N,L,tau,timeIntegration,cellmask)
tDFR3 = time.time()-t0

solutionPoints.writeSP(p,SPchoice,z12)
t0 = time.time()
x0DFR1, sol0DFR4, xDFR4, solDFR4, niterDFR4, l2DFR4 = DFR_1D.main(p,CFL,Tfin,c,D,init,grad_init,bcond,yL,yR,N,L,tau,timeIntegration,cellmask)
tDFR4 = time.time()-t0



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
plt.plot(xDFR1,solDFR1)
plt.plot(xDFR2,solDFR2,'--')
plt.plot(xDFR3,solDFR3,'-.')
plt.plot(xDFR4,solDFR4)
plt.legend(['Initial','$z_1 = z_{1 SP}$','$z_1 = z_{1 DG}$','$z_1 = $'+str(z11),'$z_1 = $'+str(z12)],loc=0)
plt.title("Solution")
plt.xlabel("x")
plt.ylabel("u")

plt.figure()
plt.plot(np.linspace(0,niterDFR1+1,niterDFR1+2)/10000.,l2DFR1)
plt.plot(np.linspace(0,niterDFR2+1,niterDFR2+2)/10000.,l2DFR2,'--')
plt.plot(np.linspace(0,niterDFR3+1,niterDFR3+2)/10000.,l2DFR3,'-.')
plt.plot(np.linspace(0,niterDFR4+1,niterDFR4+2)/10000.,l2DFR4)
plt.legend(['$z_1 = z_{1 SP}$','$z_1 = z_{1 DG}$','$z_1 = $'+str(z11),'$z_1 = $'+str(z12)],loc=0)
plt.title("L2-norm of the solution")
plt.xlabel("Iteration (10^4)")
plt.ylabel("L2-norm")


plt.show()


