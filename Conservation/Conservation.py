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

# Definition des constantes    
p = 2                                                # Degre
Ncells = 3
c = 10.
xi = solPointGen(p)                                     # SP pour le degre p
xi2 = np.insert(xi,0,-1)    
xi2 = np.append(xi2,1)                                  # SP pour le degre p+2
x = np.linspace(-1,1,501)                               # Affichage des fonctions 'continues'

xc = x[0] + (x[-1]-x[0])/Ncells*np.linspace(0,Ncells,Ncells+1)

solPointMesh = np.zeros([Ncells,len(xi)])
fluxPointMesh = np.zeros([Ncells,len(xi2)])
sol_it = np.zeros([Ncells,len(xi)])
Extrap2 = Extrap2Gen(p)
Deriv2 = D2Gen(p)
Deriv22 = D2Gen2(p)
sol_it_p2 = np.zeros([Ncells,p+3])
flux_it_p2 = np.zeros([Ncells,p+3])
derivTrueCell = np.zeros([Ncells,p+3])
derivUndersampledCell = np.zeros([Ncells,p+1])


for i in range(Ncells):
    for j in range(len(xi)):
        solPointMesh[i,j] = 0.5*(xc[i+1]+xc[i]) + 0.5*(xc[i+1]-xc[i])*xi[j]
        sol_it[i,j] = solution(solPointMesh[i,j])
    for j in range(len(xi2)):
        fluxPointMesh[i,j] = 0.5*(xc[i+1]+xc[i]) + 0.5*(xc[i+1]-xc[i])*xi2[j]


for icell in range(Ncells):
    sol_it_p2[icell,:] = np.dot(Extrap2,sol_it[icell,:])
    flux_it_p2[icell,:] = c*sol_it_p2[icell,:]

    if(c>0):
        flux_it_p2[icell,0] = c*sol_it_p2[icell-1,-1]
    elif(c<0):
        flux_it_p2[icell,-1] = c*sol_it_p2[np.mod(icell+1,Ncells),0]
        
    derivTrueCell[icell] = np.dot(Deriv22,flux_it_p2[icell,:])/(0.5*(xc[icell+1]-xc[icell]))
    derivUndersampledCell[icell] = np.dot(Deriv2,flux_it_p2[icell,:])/(0.5*(xc[icell+1]-xc[icell]))

pointsDerivTrue,derivTrue = interpolation(fluxPointMesh,derivTrueCell,p+2,100,xc[1]-xc[0],Ncells)    
pointsUSDeriv,derivUS = interpolation(solPointMesh,derivUndersampledCell,p,100,xc[1]-xc[0],Ncells)    
deriv = c*dsolution(x)

plt.figure()
#plt.plot(pointsDerivTrue,derivTrue,'-g')
plt.plot(np.reshape(fluxPointMesh,Ncells*(p+3),1),np.reshape(derivTrueCell,Ncells*(p+3),1),'og')
plt.plot(np.reshape(solPointMesh,Ncells*(p+1),1),np.reshape(derivUndersampledCell,Ncells*(p+1),1),'ob')
plt.plot(x,deriv,'--k')
plt.plot(pointsUSDeriv,derivUS,'-b')
plt.legend(('Derivative of the continuous flux','Derivative of degree p+1','True derivative'), loc='best')
plt.title("Derivatives - DFR")


#################################################################################
# FR
#################################################################################


solPointMesh = np.zeros([Ncells,len(xi)])
fluxPointMesh = np.zeros([Ncells,len(xi2)])
sol_it = np.zeros([Ncells,len(xi)])
sol_it_p2 = np.zeros([Ncells,p+3])
flux_it_p2 = np.zeros([Ncells,p+3])
g2dl = g2Derivative(xi,p)
g2dr = np.flipud(g2Derivative(xi,p))
g2dl2 = g2Derivative(xi2,p)
g2dr2 = np.flipud(g2Derivative(xi2,p))
Deriv = DGen(p)
Deriv2 = D2Gen3(p)

for i in range(Ncells):
    for j in range(len(xi)):
        solPointMesh[i,j] = 0.5*(xc[i+1]+xc[i]) + 0.5*(xc[i+1]-xc[i])*xi[j]
        sol_it[i,j] = solution(solPointMesh[i,j])
    for j in range(len(xi2)):
        fluxPointMesh[i,j] = 0.5*(xc[i+1]+xc[i]) + 0.5*(xc[i+1]-xc[i])*xi2[j]


for icell in range(Ncells):
    sol_it_p2[icell,:] = np.dot(Extrap2,sol_it[icell,:])
    flux_it_p2[icell,:] = c*sol_it_p2[icell,:]

    if(c>0):
        dflux_left = c*(sol_it_p2[icell-1,-1]-sol_it_p2[icell,0])
        dflux_right = 0.
    elif(c<0):
        dflux_left = 0.
        deflux_right = c*(sol_it_p2[np.mod(icell+1,Ncells),0]-sol_it_p2[icell,-1])
        
        
    
    derivUndersampledCell[icell] = (np.dot(Deriv,flux_it_p2[icell,1:-1]) + dflux_left*g2dl + dflux_right * g2dr)/(0.5*(xc[icell+1]-xc[icell]))
    derivTrueCell[icell] = (np.dot(Deriv2,flux_it_p2[icell,1:-1]) + dflux_left*g2dl2 + dflux_right * g2dr2)/(0.5*(xc[icell+1]-xc[icell]))
        

pointsDerivTrue,derivTrue = interpolation(fluxPointMesh,derivTrueCell,p+2,100,xc[1]-xc[0],Ncells)    
pointsUSDeriv,derivUS = interpolation(solPointMesh,derivUndersampledCell,p,100,xc[1]-xc[0],Ncells)    
deriv = c*dsolution(x)

plt.figure()
#plt.plot(pointsDerivTrue,derivTrue,'-g')
plt.plot(np.reshape(fluxPointMesh,Ncells*(p+3),1),np.reshape(derivTrueCell,Ncells*(p+3),1),'og')
plt.plot(np.reshape(solPointMesh,Ncells*(p+1),1),np.reshape(derivUndersampledCell,Ncells*(p+1),1),'ob')
plt.plot(x,deriv,'--k')
plt.plot(pointsUSDeriv,derivUS,'-b')
plt.legend(('Derivative of the continuous flux','Derivative of degree p+1','True derivative'), loc='best')
plt.title("Derivatives - FR")

