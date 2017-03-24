#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
#import sympy as sp
import scipy.optimize
import scipy.signal
import scipy.interpolate
import scipy.integrate
import math as m
import matplotlib.pyplot as plt
from polynomials import *
from gauss import *

''' authors: W. Hambli & L. Ma & R. Pile & G. Raze & A. Schmalzried'''
''' 22/02/2017 '''
''' FR_1D.py : Python library for Flux Reconstruction Method '''
''' This work follows python scripts developped by J. Vanharen and V. Joncquieres, under the supervision of G. Puigt '''



def main(p,CFL,Tfin,c,D,init,grad_init,bcond,Yl,Yr,N,L,corFun=0,timeIntegration="RK6low",cellmask="Regular"):
    ''' Order of polynomial p '''
    ''' advection celerity c, diffusion coefficient D '''
    ''' Runge-Kutta coeffcients alpha '''
    ''' Mesh composed of n isoparametric cells |-------| '''
    '''                                       -1       1 '''
    # Space discretisation for RK6
    #coeffRK=RKgamma(6)
    #alpha=RKgamma2alpha(coeffRK)
   # Space discretisation for RK6 optimized with p
    if(timeIntegration=='RK6low'):
        alpha=RKalpha6optim(p)
    elif(timeIntegration=='RK4'):
        alpha,beta = RK4()


 # Space domain
    
    # DOF
    Npoints=(p+1)*N

    #Space step

    dx = cellspacing(cellmask,L,N)

    # Cells centered about -L/2 to +L/2 AS
    x=np.zeros(N+1)
    x[0]=-0.5*N*dx[0] 
    for i in range(N):
        x[i+1]=x[i]+dx[i]


# Mesh creation with gauss points on isoparametric cells

    solPoint = solPointGen(p)
    fluxPoint = np.insert(solPoint,0,-1)
    fluxPoint = np.append(fluxPoint,1)
    solPointMesh = pointMeshGen(N,p, solPoint,dx,x)
    fluxPointMesh = pointMeshGen(N,p+2, fluxPoint,dx,x)

    dxmin = min(dx)/((p+1)**2)


    #Calcul dt 
    if(c==0. and D==0.) :
        print "Eternal frost"
        exit()
    elif c==0.: # Diffusion AS
        print "Pure Diffusion"
        dt = CFL*(dxmin)**2 /D
    elif D==0.: # Advection AS
        print "Pure Advection"
        dt= CFL*dxmin/c
    else: # Advection + Diffusion AS
        print "Advection with Diffusion"
        dtadv= CFL * dxmin / c
        dtdiff=CFL*dxmin**2 /D
        dt= min([dtadv,dtdiff])
    
    print 'dt='+ str(dt)

    
    niter = int(Tfin/dt)
# Dt for the last step to reach Tfin if Tfin =/ niter*dt
    dtfin=float(Tfin-niter*dt)
    print 'dtfin='+str(dtfin)



    
    
#Initial conditions
    
    #sol=init_triangular(solPointMesh)

    sol =np.zeros([len(solPointMesh),len(solPointMesh[0])])
    for i in range(len(solPointMesh)):
        for j in range(len(solPointMesh[0])):
            if init=='Gauss': # u0 = gaussienne 
                sol[i,j]=m.exp(-20*(solPointMesh[i,j])**2)
            elif init=='Constant': #u0 = cst
                sol[i,j] = 10.0
            elif init=='Triangle': #u0 = __/\__
                sol[i,j] = 0.0 #TODO
            elif init==1: # u0 = fonction erreur
                if 0==0:
                    sol[i,j]=(1-m.erf((solPointMesh[i,j]-c*grad_init)/2))/2 
                else:
                    sol[i,j]=(1-m.erf((solPointMesh[i,j]-c*grad_init)/(2*m.sqrt(D*grad_init))))/2
            elif init=='Sine':
                sol[i,j]=m.cos(2*2*np.pi*solPointMesh[i,j]) 
    
 
# Integration
    nG = p+1
    tG,wG = gaussRule(nG)
    intLagrange = np.zeros([len(tG),p+1])

    for i in range(len(tG)):
        for j in range(p+1):
            intLagrange[i,j]=lagrange(tG[i],solPoint,j)


#Used for the runge kutta loop
    sol0 = np.copy(sol)
#Used for comp.py
    sol00 = np.copy(sol)

# Initialisation
    sol_it = np.zeros([N,p+1])          # Values on solution points
    sol_it_ex = np.zeros([N,p+3])       # Extrapolation on solution points and interface points (discontinuous solution)
    dsol_it = np.zeros([N,p+1])         # Derivative of the discontinuous solution on solution points
    dsol_it_ex = np.zeros([N,p+3])      # Extrapolation of the derivative of the discontinuous solution on solution points and interface points
    flux_it = np.zeros([N,p+1])         # Discontinuous flux on solution points
    dflux_it_cont = np.zeros([N,p+1])   # Continuous flux derivative on solution points

    delta_flux_r = np.zeros([N])        # Left flux correction
    delta_flux_l = np.zeros([N])        # Right flux correction


    # Preparation for extrapolation (outside the loop: doesn't change inside)
    Extrap = Extrap2Gen(p)              # Lagrange extrapolation matrix on flux points
    Deriv = DGen(p)                     # Lagrange derivatives on solution points
    
    # Correction functions
    if corFun == 0.:
        dgl = g2Derivative(solPoint,p+1)            # Left correction function derivative on the solution points
        dgr = -np.flipud(dgl)                       # Right correction function derivative on the solution points
    elif corFun == 1.:
        dgl, dgr = glgrDerivatives(solPoint,p+1)    # Left and right correction function derivative on the solution points 

    if(timeIntegration=="RK4"):
        kFlux = np.zeros([4,N,p+1])                  # Stored flux at different stages

    l2 = np.zeros(niter+2)
    l2[0] = gaussIntegration(sol**2,intLagrange,wG,dx)

    ########################################################
    #                                                      #
    #                     Time Loop                        #
    #                                                      #
    ########################################################
    
    for itime in range(niter+1):

        if itime==niter:
            dti=dtfin
        else:
            dti=dt

        sol0 =np.copy(sol)
        
        for ik in range(len(alpha)):

            sol_it = np.copy(sol)


            # Extrapolation of solution on interfaces
            for icell in range(0,N):        
                sol_it_ex[icell,:] = np.dot(Extrap,sol_it[icell,:])
                if(D!=0.):
                    dsol_it[icell,:] = 2/dx[icell]*np.dot(Deriv,sol_it[icell,:])
                    dsol_it_ex[icell,:] = np.dot(Extrap,dsol_it[icell,:])

                flux_it[icell,:] = c*sol_it[icell,:] - D*dsol_it[icell,:]


            # Correction flux on interfaces computation (boundary conditions)
            if(bcond==0):       # Dirichlet BC
                if(D==0.):
                    if(c>0.):
                        delta_flux_l[0] = c*(Yl-sol_it_ex[0,0])
                        delta_flux_l[-1] = c*(sol_it_ex[-2,-1]-sol_it_ex[-1,0])
                    elif(c<0.):
                        delta_flux_r[0] = c*(sol_it_ex[1,0]-sol_it_ex[0,-1])
                        delta_flux_r[-1] = c*(Yr-sol_it_ex[-1,-1])
                else:
                    if(c>0.):
                        delta_flux_l[0] = c*(Yl-sol_it_ex[0,0])
                    delta_flux_r[0] = 0.5*( c*(sol_it_ex[1,0]-sol_it_ex[0,-1]) - D*(dsol_it_ex[1,0]-dsol_it_ex[0,-1]) )

                    delta_flux_l[-1] = 0.5*( c*(sol_it_ex[-2,-1]-sol_it_ex[-1,0]) - D*(dsol_it_ex[-2,-1]-dsol_it_ex[-1,0]) )
                    if(c<0.):
                        delta_flux_r[-1] = c*(Yr-sol_it_ex[-1,-1])

            elif(bcond==1):     # Periodic BC
                if(D==0.):
                    if(c>0.):
                        delta_flux_l[0] = c*(sol_it_ex[-1,-1]-sol_it_ex[0,0])
                        delta_flux_l[-1] = c*(sol_it_ex[-2,-1]-sol_it_ex[-1,0])
                    elif(c<0.):
                        delta_flux_r[0] = c*(sol_it_ex[1,0]-sol_it_ex[0,-1])
                        delta_flux_r[-1] = c*(sol_it_ex[0,0]-sol_it_ex[-1,-1])
                else:
                    delta_flux_l[0] = 0.5*( c*(sol_it_ex[-1,-1]-sol_it_ex[0,0]) - D*(dsol_it_ex[-1,-1]-dsol_it_ex[0,0]) )
                    delta_flux_r[0] = 0.5*( c*(sol_it_ex[1,0]-sol_it_ex[0,-1]) - D*(dsol_it_ex[1,0]-dsol_it_ex[0,-1]) )

                    delta_flux_l[-1] = 0.5*( c*(sol_it_ex[-2,-1]-sol_it_ex[-1,0]) - D*(dsol_it_ex[-2,-1]-dsol_it_ex[-1,0]) )
                    delta_flux_r[-1] = 0.5*( c*(sol_it_ex[0,0]-sol_it_ex[-1,-1]) - D*(dsol_it_ex[0,0]-dsol_it_ex[-1,-1]) )

            # Continuous flux on solution points computation (boundary cells)
            dflux_it_cont[0,:] = 2/dx[0]*(np.dot(Deriv,flux_it[0,:]) + delta_flux_l[0]*dgl + delta_flux_r[0]*dgr)
            dflux_it_cont[-1,:] = 2/dx[-1]*(np.dot(Deriv,flux_it[-1,:]) + delta_flux_l[-1]*dgl + delta_flux_r[-1]*dgr)


            # Correction flux on interfaces computation (internal cells)
            for icell in range(1,N-1):

                if(D==0.):
                    if(c>0.):
                        delta_flux_l[icell] = c*(sol_it_ex[icell-1,-1]-sol_it_ex[icell,0])
                    elif(c<0.):
                        delta_flux_r[icell] = c*(sol_it_ex[np.mod(icell+1,N),0]-sol_it_ex[icell,-1])
                else:
                    delta_flux_l[icell] = 0.5*( c*(sol_it_ex[icell-1,-1]-sol_it_ex[icell,0]) - D*(dsol_it_ex[icell-1,-1]-dsol_it_ex[icell,0]) )
                    delta_flux_r[icell] = 0.5*( c*(sol_it_ex[np.mod(icell+1,N),0]-sol_it_ex[icell,-1]) - D*(dsol_it_ex[np.mod(icell+1,N),0]-dsol_it_ex[icell,-1]) )
                
                # Continuous flux on solution points computation (internal cells)
                dflux_it_cont[icell,:] = 2/dx[icell]*(np.dot(Deriv,flux_it[icell,:]) + delta_flux_l[icell]*dgl + delta_flux_r[icell]*dgr)



            # Solution update
            sol = sol0 - dti * alpha[ik]*dflux_it_cont
            if(timeIntegration=="RK4"):         # Saving stages
                kFlux[ik,:,:] = dflux_it_cont
        
        l2[itime+1] = gaussIntegration(sol**2,intLagrange,wG,dx)

        if(timeIntegration=="RK4"):             # Combining stages
            sol = sol0
            for ik in range(len(beta)):
                sol = sol - dti*beta[ik]*kFlux[ik,:,:]

        if itime==niter:
            dt1=dt
            print "-----------------------------------------------"
            print "Iteration:"+str(itime) +",Time: " + str(itime*dt1+dtfin) + "/" + str(niter*dt+dtfin)
        else:
            if divmod(itime,1000)[1]==0:
                print "-----------------------------------------------"
                print "Iteration:"+str(itime)+ ",Time: " + str((itime + 1)*dt) + "/" + str(niter*dt+dtfin)

    print "-----------------------------------------------"

 # to reshape the matrix into a vector '''

    solPointMesh = solPointMesh.reshape((p + 1) * N)
    sol = sol.reshape(((p + 1) * N))
    sol00 = np.copy(sol)
    solPointMesh00=np.copy(solPointMesh)
    
# final number of points for interpolation
    h=1000
    solPointMesh,sol=interpolation(solPointMesh,sol,p,h,dx[0],N)
    # solPointMesh00,sol00=interpolation(solPointMesh00,sol00,p,h,dx[0],N)        #-> Il y a l'air d'avoir un souci avec l'interpolation, ca ne donne pas la bonne solution de depart (GR)
    
    

    return solPointMesh00, sol00, solPointMesh, sol, niter, l2



'''Part 1 Position of points and mesh of the domain'''

def solPointGen(p):
    ''' Compute solution points for an isoparametric cell with p + 1 Gauss-Lobatto points '''
    solPoint = np.zeros(p+1)
    for i in range(len(solPoint)):
        solPoint[i] = - np.cos(np.pi * (2. * (i + 1) - 1) / (2 * (p + 1))) # Peut-être qu'il faut faire solPoint[i+1], not sure --> Nope, c'est correct car en python on commence à 0
    return solPoint 

    
def pointMeshGen(N,p, point,dx,xreal):
    ''' Compute flux or solution points '''
    ''' for the whole mesh composed of n cells '''
    ''' in real domain '''

    Jac=dx/2
    x = np.zeros([N,p+1])
    for i in range(N):
        for j in range(p+1):
            x[i,j] = point[j]*Jac[i] +(xreal[i]+xreal[i+1])/2 
           

    return x


def init_triangular(x):
    y=np.zeros([len(x),len(x[0])])
    for i in range(len(x)):
        for j in range(len(x[0])):
            if (divmod(np.floor(x[i,j]),2)[1]==0):
                y[i,j]=x[i,j]-np.floor(x[i,j])
            else:
                y[i,j]=-x[i,j]+np.floor(x[i,j])+1
    return y


''' Part 2 extrapolation  of the flux F'''


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
            Extrap[i,j]=lagrange(solPoint[i],solPoint,j) # row = lagrange polynomial value on solpoint i | column = lagrange polynomial non zero on j
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
            Extrap[i,j]=lagrange(fluxPoint[i],solPoint,j) # row = LagPol value on solpoint i | column = LagPol from j
    #print Extrap[:,:]
    return Extrap


''' Part 3 Riemann solver : Godunov '''

def f(u,c):
    res=c*u
    return res
    

def fprime(u,c):
    res=c;
    return res;

def Godunov(uL,uR,c):
    res=0.0
    if(uL==uR):
        res=f(uL,c);
    else:
        if (fprime(uL,c)>fprime(uR,c)):
            sigma=(f(uR,c)-f(uL,c))/(uR-uL)
            if(sigma==0):
                res=f(uL,c);
            else:
                if(sigma>0.0):
                    res=f(uL,c);
                else:
                    res=f(uR,c);
        else:
            if(fprime(uL,c)>0.0):
                res=f(uL,c);
            else:
                if(fprime(uR,c)<0.0):
                    res=f(uR,c);
                else:
                    res=0.0;
    return res;


''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

''' Part 4 Derivation of the flux'''
def DGen(p, phi=1.0):    #RP
    ''' Compute d/dx discretization with '''
    ''' the Spectral Difference method order p '''
    ''' phi = 1.0 <=> upwind  flux '''
    ''' phi = 0.0 <=> centred flux '''

    solPoint = solPointGen(p)

    # Addition of interface points 

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
    fluxPoint = np.insert(solPoint,0,-1) 
    fluxPoint = np.append(fluxPoint,1)

    ''' Compute the derivatives matrix '''
    D = np.zeros([p+3, p+1])
    for i in range(p+3):
        for j in range(p+1):
            D[i, j] = lagrangeDerivative(fluxPoint[i], j, solPoint) # row = d(LagPol)/dx value on solpoint i | column = d(LagPol from j)/dx
    return D


''' Part 5 Runge Kutta'''

def RKgamma2alpha(gamma):
    ''' Transformation from 'gamma' to 'alpha' Runge-Kutta coefficients '''
    alpha = np.zeros(len(gamma))
    prod = 1.
    for i in range(len(gamma)):
        alpha[-i-1] = gamma[i]/prod
        prod = prod*alpha[-i-1]
    return alpha


def RKgamma(order):
    ''' Runge-Kutta coefficients for time integration '''
    gamma = np.zeros(order)
    for i in range(order):
        gamma[i] = 1. / np.math.factorial(i + 1)
    return gamma

def RKalpha6optim(p):
    ''' Runge-Kutta coefficients for time integration optimized for order 6'''
    alpha = np.zeros(6)
    alpha[2]=0.24662360430959
    alpha[3]=0.33183954253762
    alpha[4]=0.5
    alpha[5]=1.0
    if (p==2):
            alpha[0]=0.05114987425612
            alpha[1]=0.13834878188543
    if (p==3):
            alpha[0]=0.07868681448952
            alpha[1]=0.12948018884941
    if (p==4):
            alpha[0]=0.06377275785911
            alpha[1]=0.15384606858263
    if (p==5):
            alpha[0]=0.06964990321063
            alpha[1]=0.13259436863348
    if (p==6):
            alpha[0]=0.06809977676724
            alpha[1]=0.15779153065865
    if (p==7):
            alpha[0]=0.06961281995158
            alpha[1]=0.14018408222804
    if (p==8):
            alpha[0]=0.07150767268798
            alpha[1]=0.16219675431011
    if (p==9):
            alpha[0]= 0.06599710352324
            alpha[1]=0.13834850670675
    if (p==10):
            alpha[0]=0.07268810031422
            alpha[1]=0.16368178688643
    return alpha

def RK4():
    alpha = np.zeros(4)
    alpha[0] = 0.5 
    alpha[1] = 0.5
    alpha[2] = 1.

    beta = np.zeros(4)
    beta[0] = 1.0/6.0
    beta[1] = 1.0/3.0
    beta[2] = 1.0/3.0
    beta[3] = 1.0/6.0

    return alpha, beta

''' Part 6 Post Processing'''
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
    solPointMesh = solPointMesh.reshape((int(len(solPointMesh)/(p+1)),int(p+1)))
    sol = sol.reshape((int(len(sol)/(p+1)),int(p+1)))    
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
    
def display(x,y,p,dx,N):
    '''Display y(x) with a higher resolution'''
    h=1000
    xx,yy=interpolation(x,y,p+2,h,dx,N)
    plt.plot(xx,yy, 'r-')

def cellspacing(maskoption,L,N):
    mask = np.ones(N)   #irregular cell size mask proportional to regular size
    dx_reg = L/N # Regular cell spacing (AS)
    if(maskoption=='Irregular'):
        for i in range(N): #mask customization
            if(i%2==0):
                mask[i] = mask[i]*1.1
            else:
                mask[i] = mask[i]*0.9
            
    dx = mask*dx_reg
    dx = dx*L/np.sum(dx) #normalisation to match the length of the domain
    return dx
    
