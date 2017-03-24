# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 18:10:46 2017

@author: Ghislain
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 14:38:58 2017

@author: Ghislain
"""

import numpy as np
import matplotlib.pyplot as plt
from polynomials import *
from gauss import *

def solPointGen(p):
    ''' Compute solution points for an isoparametric cell with p + 1 Gauss-Lobatto points '''
    solPoint = np.zeros(p+1)
    for i in range(len(solPoint)):
        solPoint[i] = - np.cos(np.pi * (2. * (i + 1) - 1) / (2 * (p + 1))) # Peut-être qu'il faut faire solPoint[i+1], not sure --> Nope, c'est correct car en python on commence à 0
    return solPoint 


def lagrange(x,xi,i):
    ''' Lagrange polynomial for extrapolation 
        zeros on xi except for i '''
    res = 1.0
    for s in range(len(xi)):
        if i != s:
            res = res * (x - xi[s]) / (xi[i] - xi[s])
    return res


def lagrangeDerivative(x, i, xi):
    ''' Lagrange polynomial derivative at x 
        xi : zeros of polynomial
        i = unimodular non-zero exception point index '''
    res = 1.0
    somme = 0.0
    if x in xi:
        ind = np.linspace(0,len(xi)-1,len(xi))
        ind = int(ind[x==xi])
        for s in range(len(xi)):
            if i != s and ind != s:
                res = res * (x - xi[s]) / (xi[i] - xi[s])
                somme = somme + 1.0/(x -xi[s])
            elif i != s:
                res = res/ (xi[i] - xi[s])
        if i != ind:        
            der = res*(x-xi[ind])*somme + res
        else:
            der = res*somme
                
    else:
        for s in range(len(xi)):
            if i != s:
                res = res * (x - xi[s]) / (xi[i] - xi[s])
                somme = somme + 1.0/(x -xi[s])
        der = somme*res
        
    return der
 

def evaluateLagrange(x,xi):

    lag = np.zeros([len(x),len(xi)])
    dlag = np.zeros([len(x),len(xi)])
    
    for j in range(len(x)):
        for i in range(len(xi)):
            lag[j][i] = lagrange(x[j],xi,i)
            dlag[j][i] = lagrangeDerivative(x[j],i,xi)
            
    return lag, dlag

def ExtrapGen(p): # Fonction sans intérêt : on extrapole sur les points solutions -> matrice diagonale !
    ''' Extrapolation through lagrange polynomials on p+1 points '''
    ns = p+1
    solPoint = solPointGen(p) # A optimiser : on peut le generer 1! fois AS

    #RP nf = p+2
    #RP fluxPoint = fluxPointGen(p)

    Extrap=np.zeros([ns,ns]); #(RP) ns + 1 ? -> on extrapole seulement apd des p+1 pour l'instant... AS
    for i in range(ns):
        for j in range(ns):
            #Extrap[i,j]=lagrange(fluxPoint[i],solPoint,j);
            Extrap[i,j]=lagrange(solPoint[i],solPoint,j); # row = lagrange polynomial value on solpoint i | column = lagrange polynomial non zero on j
    return Extrap
 
def Extrap2Gen(p): 
    ''' Extrapolation matrix for the P+3 reconstruction with border terms ''' 
    solPoint = solPointGen(p) 

    #Ajout des bords aux points solutions
    fluxPoint = np.insert(solPoint,0,-1) 
    fluxPoint = np.append(fluxPoint,1)

    Extrap=np.zeros([p+3,p+1]); #Extrapolation of p+1 LagPol to p+3 points
    for i in range(p+3):
        for j in range(p+1):
            Extrap[i,j]=lagrange(fluxPoint[i],solPoint,j); # row = LagPol value on solpoint i | column = LagPol from j
    #print Extrap[:,:]
    return Extrap
    

def DGen(p, phi=1.0):    #RP
    ''' Compute d/dx discretization with '''
    ''' the Spectral Difference method order p '''
    ''' phi = 1.0 <=> upwind  flux '''
    ''' phi = 0.0 <=> centred flux '''

    solPoint = solPointGen(p)

    ''' Compute the derivatives matrix '''
    D = np.zeros([p+1, p+1])
    for i in range(p+1):
        for j in range(p+1):
            D[i, j] = lagrangeDerivative(solPoint[i], j, solPoint) # row = d(LagPol)/dx value on solpoint i | column = d(LagPol from j)/dx
    return D    

def D2Gen(p, phi=1.0):    #RP
    ''' Compute d/dx discretization with '''
    ''' the Spectral Difference method order p '''
    ''' phi = 1.0 <=> upwind  flux '''
    ''' phi = 0.0 <=> centred flux '''

    solPoint = solPointGen(p)

    # Addition of interface points 
    solPoint = np.insert(solPoint,0,-1) 
    solPoint = np.append(solPoint,1)

    ''' Compute the derivatives matrix '''
    D = np.zeros([p+1, p+3])
    for i in range(p+1):
        for j in range(p+3):
            D[i, j] = lagrangeDerivative(solPoint[i+1], j, solPoint) # row = d(LagPol)/dx value on solpoint i | column = d(LagPol from j)/dx
    return D
    
def D2Gen2(p, phi=1.0):    #RP
    ''' Compute d/dx discretization with '''
    ''' the Spectral Difference method order p '''
    ''' phi = 1.0 <=> upwind  flux '''
    ''' phi = 0.0 <=> centred flux '''

    solPoint = solPointGen(p)

    # Addition of interface points 
    solPoint = np.insert(solPoint,0,-1) 
    solPoint = np.append(solPoint,1)

    ''' Compute the derivatives matrix '''
    D = np.zeros([p+3, p+3])
    for i in range(p+3):
        for j in range(p+3):
            D[i, j] = lagrangeDerivative(solPoint[i], j, solPoint) # row = d(LagPol)/dx value on solpoint i | column = d(LagPol from j)/dx
    return D    

def D2Gen3(p, phi=1.0):    #RP
    ''' Compute d/dx discretization with '''
    ''' the Spectral Difference method order p '''
    ''' phi = 1.0 <=> upwind  flux '''
    ''' phi = 0.0 <=> centred flux '''

    solPoint = solPointGen(p)
    # Addition of interface points 
    solPoint = np.insert(solPoint,0,-1) 
    solPoint = np.append(solPoint,1)

    ''' Compute the derivatives matrix '''
    D = np.zeros([p+3, p+1])
    for i in range(p+3):
        for j in range(p+1):
            D[i, j] = lagrangeDerivative(solPoint[i], j, solPoint[1:-1]) # row = d(LagPol)/dx value on solpoint i | column = d(LagPol from j)/dx
    return D  

def InterpGen(p,h,xinterp,solPointMesh):

    ns = len(solPointMesh[0])
   

    nf = len(xinterp[0])
    

    Interp=np.zeros([nf,ns]);
    for i in range(nf):
        for j in range(ns):
            Interp[i,j]=lagrange(xinterp[0,i],solPointMesh[0,:],j);
   
    return Interp

def interpolation(solPointMesh, sol, p, h,dx,Ncell):
    Extrap=ExtrapGen(p)
    
    J=dx/2
    #solPointMesh = solPointMesh.reshape((int(len(solPointMesh)/(p+1)),int(p+1)))
    #sol = sol.reshape((int(len(sol)/(p+1)),int(p+1)))    
    xinterp = np.zeros((solPointMesh.shape[0], len(np.arange(0, dx+0.000000000001, dx/h))))
    for icell in range(len(solPointMesh)):
        a=np.arange(dx*icell, dx*(icell+1)+0.00000000001, dx/h)-dx*Ncell/2
        xinterp[icell,:] = a
        
    
    Interp=InterpGen(p,h,xinterp,solPointMesh)
    yinterp = 0.*xinterp
    
    
    
    for icell in range(len(solPointMesh)):
        yinterp[icell,:]=np.dot(Interp,sol[icell,:])


    
    xinterp = xinterp.reshape(xinterp.shape[0]*xinterp.shape[1])
    yinterp = yinterp.reshape(yinterp.shape[0]*yinterp.shape[1])
    return xinterp, yinterp
    
   
def solution(x):
    return np.exp(-(x)**2/0.1)
    

def dsolution(x):
    return -2*x/0.1*np.exp(-(x)**2/0.1)    


SP = np.array([0.5])
p = len(SP)-1
nG = np.ceil(0.5*(p+1))
wG,tG = gaussRule(nG)

FP = np.append(SP,1)
FP = np.insert(FP,0,-1)
yi = np.array([2.,1.,3.])

print FP
print yi
