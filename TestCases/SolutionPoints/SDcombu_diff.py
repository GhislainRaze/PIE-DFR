import numpy as np
import sympy as sp
import scipy.optimize
import scipy.signal
import scipy.interpolate
import scipy.integrate
import math as m
import matplotlib.pyplot as plt
from gauss import *
''' author: V. Joncquieres '''
''' 23/03/15 '''
''' SD combu.py : Python library for Spectral Difference Method in combustion '''

def main(p,method,CFL,Tfin,c,D,init,grad_init,bcondL,bcondR,Yl,Yr,N,L,timeIntegration="RK6low",cellmask="Regular"):
    ''' Order of polynomial p '''
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
    
    dx1 = np.min(dx)/(p+1)
#Calcul dt 
    if (method==0):
#test in case
        if c==0.:
            print "!!! c must have a value !!!"
            exit()
        else:
            dt = CFL * dx1 / (c*(p+1))
    elif(method==1):
        if D==0.:
            print "!!! D must have a value !!!"
            exit()
        else:
            dt= CFL*(dx1/(p+1))**2 /D
    else:
        if (D==0.) or (c==0.):
            print "!!!c and D must have  values !!!"
            exit()
        else:
            dtadv= CFL * min(dx) / (c*(p+1))
            dtdiff=CFL*(min(dx)/(p+1))**2 /D
            dt= min([dtadv,dtdiff])
    
    print 'dt='+ str(dt)
    
    niter = int(Tfin/dt)
# Dt for the last step to reach Tfin if Tfin =/ niter*dt
    dtfin=float(Tfin-niter*dt)
    print 'dtfin='+str(dtfin)


#Mesh creation with gauss points on isoparametric cells

    solPoint = solPointGen(p)
    solPointMesh = pointMeshGen(N,p, solPoint,dx,x)
#Initialisation
    

    #sol=init_triangular(solPointMesh)

    sol =np.zeros([len(solPointMesh),len(solPointMesh[0])])
    for i in range(len(solPointMesh)):
        for j in range(len(solPointMesh[0])):
            if init==0:
                sol[i,j]=m.exp(-20*(solPointMesh[i,j])**2)
            elif init==1: 
                if method==0:
                    sol[i,j]=(1-m.erf((solPointMesh[i,j]-c*grad_init)/2))/2
                else:
                    sol[i,j]=(1-m.erf((solPointMesh[i,j]-c*grad_init)/(2*m.sqrt(D*grad_init))))/2
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

    sol_flux_it=np.zeros([N,p+2])
    grad=np.zeros([N,p+1])
    gradF=np.zeros([N,p+2])
    lapl=np.zeros([N,p+1])
    Dadv= np.zeros([N, p + 1])


# preparation for extrapolation (outside the loop: doesn't change inside)

    
    Extrap = ExtrapGen(p)
    Deriv=DGen(p)

    if(timeIntegration=="RK4"):
        kFlux = np.zeros([4,N,p+1])                  # Stored flux at different stages


    l2 = np.zeros(niter+2)
    l2[0] = gaussIntegration(sol**2,intLagrange,wG,dx)

#Time loop until niter+1 for the last iteration to reach Tfin
    for itime in range(niter+1):
       
        if itime==niter:
            dti=dtfin
        else:
            dti=dt


        sol0 =np.copy(sol)
#Range Kutta loop '''
        for ik in range(len(alpha)):
#loop on number of cells'''
            for icell in range(N):
                
# real domain to isoparametric dimension 
                #Jacobian in 1D
                J=dx[0]/2
                
                
                sol_it= sol/J
# Extrapolation of solution on fluxes points
                
                sol_flux_it[icell,:]=np.dot(Extrap,sol_it[icell,:])
                if bcondL==0:
                    sol_flux_it[0,0]=Yl/J
                elif bcondL==1:
                    sol_flux_it[0,0]=sol_flux_it[-1,-1]	
                if bcondR==0:
                    sol_flux_it[-1,-1]=Yr/J
                elif bcondL==1:
                    sol_flux_it[-1,-1]=sol_flux_it[0,0]
		# ADVECTION  Computation of fluxes at flux points : for advection F=cY
                
                if (method!=1):
                    flux_it=c*sol_flux_it

# isoparametric dimension to real domain not necessary in 1D
                    flux = flux_it
#*J**(-1)*J

# Treatment of interfaces with upwind scheme  for Advection
                    for k in range(1,N):
                        flux[k,0]=flux[k-1,-1]
               
                
                

# real domain to isoparametric transformation non necessary in 1D
                #flux_it = flux
#/(J**(-1)*J)
                
#Derivation of sol on flux point with Deriv 
   
                    Dadv[icell, :] = np.dot(Deriv, flux_it[icell,:])

                #DIFFUSION

# Treatment of interfaces with mean for diffusion
                if (method !=0):
                    sol_flux_it_tmp = np.copy(sol_flux_it)
                    for k in range(1,N):
                        sol_flux_it[k,0] = 0.5*(sol_flux_it_tmp[k,0]+sol_flux_it_tmp[k-1,-1]) 
                    for k in range(0,N-1):
                        sol_flux_it[k,-1]= 0.5*(sol_flux_it_tmp[k+1,0]+ sol_flux_it_tmp[k,-1])
                
# Derivation of sol on flux point give a grad on sol point 

                    grad[icell, :] = np.dot(Deriv, sol_flux_it[icell,:])
                    gradF[icell,:] = np.dot(Extrap,grad[icell,:])/J
                
# Treatment of interfaces with mean for diffusion
                    grad_tmp = np.copy(gradF)
                    for k in range(1,N):
                        gradF[k,0] = 0.5*(grad_tmp[k,0]+grad_tmp[k-1,-1]) 
                    for k in range(0,N-1):
                        gradF[k,-1]= 0.5*(grad_tmp[k+1,0]+ grad_tmp[k,-1])
                

# Derivation of grad
                    lapl[icell, :] = np.dot(Deriv, gradF[icell,:])
            
            # Correction: for periodic BCs (GR)
            Dadv[0, :] = np.dot(Deriv, flux_it[0,:])
            lapl[0, :] = np.dot(Deriv, gradF[0,:])   

            sol = sol0 + dti *  alpha[ik]*(D*lapl-Dadv)
            if(timeIntegration=="RK4"):         # Saving stages
                kFlux[ik,:,:] = Dadv-D*lapl

        if(timeIntegration=="RK4"):             # Combining stages
            sol = sol0
            for ik in range(len(beta)):
                sol = sol - dti*beta[ik]*kFlux[ik,:,:]
                
        l2[itime+1] = gaussIntegration(sol**2,intLagrange,wG,dx)
        if itime==niter:
            dt1=dt
            print "Iteration:"+str(itime) +",Time: " + str(itime*dt1+dtfin) + "/" + str(niter*dt+dtfin)
        else:
            if divmod(itime,1000)[1]==0:
                print "-----------------------------------------------"
                print "Iteration:"+str(itime)+ ",Time: " + str((itime + 1)*dt) + "/" + str(niter*dt+dtfin)

    




# to reshape the matrix into a vector '''

    solPointMesh = solPointMesh.reshape((p + 1) * N)
    sol = sol.reshape(((p + 1) * N))
    sol00 = sol00.reshape((p + 1) * N)
    solPointMesh00=np.copy(solPointMesh)
    
# final number of points for interpolation
    h=1000
    solPointMesh,sol=interpolation(solPointMesh00,sol,p,h,dx[0],N)
    solPointMesh00,sol00=interpolation(solPointMesh00,sol00,p,h,dx[0],N)
    return solPointMesh00, sol00,solPointMesh, sol, niter, l2





''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''



'''Part 1 Position of points and mesh of the domain'''

def solPointGen(p):
    return np.loadtxt("SP.txt")


def pointMeshGen(N,p, point,dx,xreel):
    ''' Compute flux or solution points '''
    ''' for the whole mesh composed of n cells '''
    ''' in real domain '''

    Jac=dx/2
    x = np.zeros([N,p+1])
    for i in range(N):
        for j in range(p+1):
            x[i,j] = point[j]*Jac[i] +(xreel[i]+xreel[i+1])/2 
           

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

def fluxPointGen(p):
    ''' Compute flux points : -1, 1 and roots of Legendre polynomials for
    an isoparametric cell '''
    fluxPoint = np.zeros(p+2)
    fluxPoint[0] = -1.0
    fluxPoint[-1] = 1.0
    c = np.zeros(p + 1)
    c[-1] = 1.0
    coeff = np.polynomial.legendre.legroots(c)
    for i in range(len(coeff)):
        fluxPoint[1 + i] = coeff[i]
    return fluxPoint

def lagrange(x,xi,i):
    ''' Lagrange polynomial for extrapolation'''
    res = 1.0
    for s in range(len(xi)):
        if i != s:
            res = res * (x - xi[s]) / (xi[i] - xi[s])
    return res

def ExtrapGen(p):

    ns = p+1
    solPoint = solPointGen(p)

    nf = p+2
    fluxPoint = fluxPointGen(p)

    Extrap=np.zeros([nf,ns]);
    for i in range(nf):
        for j in range(ns):
            Extrap[i,j]=lagrange(fluxPoint[i],solPoint,j);
   
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



def lagrangeDerivative(x, i, xi):
    ''' Lagrange polynomial derivative '''
    res = 1.0
    xvar = sp.symbols("xvar")
    for s in range(len(xi)):
        if i != s:
            res = res * (xvar - xi[s]) / (xi[i] - xi[s])
    dev = sp.diff(res, xvar)
    return dev.subs(xvar, x)

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

''' Part 4 Derivation of the flux'''

def DGen(p, phi=1.0):
    ''' Compute d/dx discretization with '''
    ''' the Spectral Difference method order p '''
    ''' phi = 1.0 <=> upwind  flux '''
    ''' phi = 0.0 <=> centred flux '''

    ns = p+1
    solPoint = solPointGen(p)

    nf = p+2
    fluxPoint = fluxPointGen(p)

    ''' Compute the derivatives matrix '''
    D = np.zeros([ns, nf])
    for i in range(ns):
        for ii in range(nf):
            D[i, ii] = lagrangeDerivative(solPoint[i], ii, fluxPoint)


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
