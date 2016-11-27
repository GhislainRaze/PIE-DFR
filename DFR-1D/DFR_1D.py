import numpy as np
#import sympy as sp
import scipy.optimize
import scipy.signal
import scipy.interpolate
import scipy.integrate
import math as m
import matplotlib.pyplot as plt
''' author: R. Pile & A. Schmalzried '''
''' 19/11/2016 '''
''' DFR_1D.py : Python library for direct Flux reconstruction Method '''



def main(p,method,CFL,Tfin,c,D,init,grad_init,bcondL,bcondR,Yl,Yr):
    ''' Order of polynomial p '''
    ''' advection celerity c, diffusion coefficient D '''
    ''' Runge-Kutta coeffcients alpha '''
    ''' Mesh composed of n isoparametric cells |-------| '''
    '''                                       -1       1 '''
    # Space discretisation for RK6
    #coeffRK=RKgamma(6)
    #alpha=RKgamma2alpha(coeffRK)
   # Space discretisation for RK6 optimized with p
    alpha=RKalpha6optim(p)

 # Space domain

    # DOF
    Npoints=300
    #Number of cells
    N=int(Npoints/(p+1)+1) #In case not enough points for order p, decrease in the number of cells ? AS
    #Space step
    dx1=0.001
    dx=np.linspace(dx1*(p+1),dx1*(p+1),N) #Cell step AS

    # Cells centered about -L/2 to +L/2 AS
    L=np.sum(dx)
    x=np.zeros(N+1)
    x[0]=-N*dx[0]/2.0 
    for i in range(N):
        x[i+1]=x[i]+dx[i]
    
#Calcul dt 
    if (method==0): # Advection AS
#test in case
        if c==0.:
            print "!!! c must have a value !!!"
            exit()
        else:
            dt = CFL * (dx1/p+1) / c
    elif(method==1): # Diffusion AS
        if D==0.:
            print "!!! D must have a value !!!"
            exit()
        else:
            dt= CFL*(dx1/(p+1))**2 /D
    else: # Advection + Diffusion AS
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
    
    
#Initial conditions
    
    #sol=init_triangular(solPointMesh)

    sol =np.zeros([len(solPointMesh),len(solPointMesh[0])])
    for i in range(len(solPointMesh)):
        for j in range(len(solPointMesh[0])):
            if init==0: # u0 = gaussienne 
                sol[i,j]=m.exp(-(solPointMesh[i,j]+0.10)**2/0.0001) #Calcul des points flux. USELESS FOR DFR (RP) --> il s'agit plutot de la solution initiale (Gaussienne) ! 
            elif init==1: # u0 = fonction erreur
                if method==0:
                    sol[i,j]=(1-m.erf((solPointMesh[i,j]-c*grad_init)/2))/2 
                else:
                    sol[i,j]=(1-m.erf((solPointMesh[i,j]-c*grad_init)/(2*m.sqrt(D*grad_init))))/2
    
    #sol = solPointMesh #(RP) --> why ? solPointMesh est la matrice des noeuds en abscisse !
    
#Used for the runge kutta loop
    sol0 = np.copy(sol)
#Used for comp.py
    sol00 = np.copy(sol)


# Initialisation

#RP    sol_flux_it=np.zeros([N,p+2])
    sol_flux_it=np.zeros([N,p+1])
    sol_it = np.zeros([N,p+1])
    sol_it_p2 = np.zeros([N,p+1])
# RP   grad=np.zeros([N,p+1])
#    gradF=np.zeros([N,p+2])
#    lapl=np.zeros([N,p+1])
#    Dadv= np.zeros([N, p + 1])
    sol_flux_p2 = np.zeros([N,p+3])
    aux_var2_it = np.zeros([N,p+3])
    aux_var_it = np.zeros([N,p+3])
    flux_Dd = np.zeros([N,p+1])
    flux_Dd2 = np.zeros([N,p+3])
    flux = np.zeros([N,p+1])


# preparation for extrapolation (outside the loop: doesn't change inside)

    
    Extrap = ExtrapGen(p) # Lagrange extrapolation matrix of order P
    
    Extrap2 = Extrap2Gen(p) # Lagrange extrapolation matrix of order P+2
    
    Deriv2=D2Gen(p) # Lagrange derivative extrapolation matrix of order P+2

#Time loop until niter+1 for the last iteration to reach Tfin
    for itime in range(niter+1): # ??? On fait niter+1 itérations au lieu de niter ? AS

        if itime==niter:
            dti=dtfin
        else:
            dti=dt


        sol0 =np.copy(sol)
        
#Range Kutta loop '''
        for ik in range(len(alpha)):

            for icell in range(N): # Extrapolation loop on cell boundaries
                sol_it_p2[icell,:] = np.dot(Extrap2,sol_it[icell,:])
                sol_it_p2[0,0] = bcondL
                sol_it_p2[-1,-1] = bcondR

################################################################################################################### Version DFR ci-dessous (RP)

        # Extrapolation of solution on borders (RP)
            sol_it_tmp = np.copy(sol_it)
            
            sol_it[icell,0] = 0.5*(sol_it_tmp[icell,0]+sol_it_tmp[icell-1,p+2]) 
            sol_it[icell,p+2]= 0.5*(sol_it_tmp[icell+1,0]+ sol_it_tmp[icell,p+2])
                    
    # Lagrangian interpolation of the complete solution (RP)                
            sol_flux_p2[icell,:]=np.dot(Extrap2,sol_flux_it[icell,:])
            
    # Auxiliary variable computation (RP)
            aux_var2_it[icell,:] = np.dot(Deriv2,sol_flux_it[icell,:])
            
    # Undersampling of the auxiliary variable (RP)
            del aux_var2_it[icell,0]
            del aux_var2_it[icell,p+2]
            aux_var_it[icell,:] = np.dot(Extrap,aux_var2_it[icell,:])
            
            
    # Computation of the flux (RP) 
            #RP TODO : A function computing the flux for more complex problem
            flux_d[icell,:] = c*sol_flux_it_tmp[icell,:] + D*aux_var_it[icell,:] #RP A REMPLACER
            flux_Dd[icell,:]=np.dot(Extrap,flux_d[icell,:]) #Ligne inutile dans ce cas, mais utile si on a un pblm non linéaire
            
    # Extrapolation of the flux on borders (RP)        
            flux_Dd_temp = np.copy(flux_Dd)
            flux_Dd[icell,0] = 0.5*(flux_Dd_tmp[icell,0]+flux_Dd_tmp[icell-1,p+2]) + tau*(sol_flux_it_tmp[icell,0] - sol_flux_it_tmp[icell-1,p+2]) 
            flux_Dd[icell,p+2]= 0.5*(flux_Dd_tmp[icell+1,0]+ flux_Dd_tmp[icell,p+2]) + tau*(sol_flux_it_tmp[icell+1,0]+ sol_flux_it_tmp[icell,p+2])
            
    # Lagrangian interpolation of the complete flux (RP)
            flux_Dd2[icell,:]=np.dot(Extrap2,flux_Dd[icell,:])
            
     # Undersampling of the  flux (RP)
            del flux_Dd2[icell,0]
            del flux_Dd2[icell,p+2]
            flux[icell,:] = np.dot(Extrap,aux_var2_it[icell,:])
            
    
# TOFINISH (RP)


################################################################################################################### Version SD ci-dessous (RP)

''' bullship spoiler down below !
#loop on number of cells
            for icell in range(N):
                
# real domain to isoparametric dimension 
                #Jacobian in 1D
                J=dx[0]/2 # pussy -> c'est dx[icell] qu'il faut utiliser !
                sol_it= sol/J #cela n'a pas d'interet ! c'est sur les distances qu'il faut utiliser le jacobien...
                
                
# Extrapolation of solution on fluxes [= solutions (RP)] points -> ne sert à rien ! car les points sont confondus avec les solpoint, donc valent zéro partout sauf au point considéré
                sol_flux_it[icell,:]=np.dot(Extrap,sol_it[icell,:]) # Si c'est bien le produit matriciel que je pense alors les 
                                                                    #lagrangiens doivent etre ranges de manière à avoir les valeurs 
                                                                    #en les points intérieurs sur les colonnes de Extrap cf en dessous AS
                
    # Les CL sont à l'interieur de la boucle parce que sinon on a un problème de dimensions avec la ligne ci-dessus (RP)
                if bcondL==0:
                    sol_flux_it[0,0]=Yl/J
                elif bcondL==1:
                    sol_flux_it[0,0]=sol_flux_it[-1,-1]	
                if bcondR==0:
                    sol_flux_it[-1,-1]=Yr/J
                elif bcondL==1:
                    sol_flux_it[-1,-1]=sol_flux_it[0,0]
                    
   
    # MAYBE USELESS : TOCHECK (RP)                  
#		# ADVECTION  Computation of fluxes at flux points : for advection F=cY
#                
#                if (method!=1):
#                    flux_it=c*sol_flux_it
#
## isoparametric dimension to real domain not necessary in 1D
#                    flux = flux_it
##*J**(-1)*J
'''


# Treatment of interfaces with upwind scheme  for Advection
                    for k in range(1,N):
                        flux[k,0]=flux[k-1,-1] # Continuité du flux à gauche. Probablement pas utile pour DFR (RP)
               
                
# real domain to isoparametric transformation non necessary in 1D
                #flux_it = flux
#/(J**(-1)*J)


#Derivation of sol on flux point with Deriv 
   
                    Dadv[icell, :] = np.dot(Deriv, flux_it[icell,:]) # calcul de df/dx , soit notre q

                #DIFFUSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
                
            sol = sol0 + dti *  alpha[ik]*(D*lapl-Dadv)
        
                
        
        if itime==niter:
            dt1=dt
            print "Iteration:"+str(itime) +",Time: " + str(itime*dt1+dtfin) + "/" + str(niter*dt+dtfin)
        else:
            if divmod(itime,10)[1]==0:
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
    return solPointMesh00, sol00,solPointMesh, sol, niter





''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''



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

#def fluxPointGen(p):
#    ''' Compute flux points : -1, 1 and roots of Legendre polynomials for
#    an isoparametric cell '''
#    fluxPoint = np.zeros(p+2)
#    fluxPoint[0] = -1.0
#    fluxPoint[-1] = 1.0
#    c = np.zeros(p + 1)
#    c[-1] = 1.0
#    coeff = np.polynomial.legendre.legroots(c)
#    for i in range(len(coeff)):
#        fluxPoint[1 + i] = coeff[i]
#    return fluxPoint

def lagrange(x,xi,i):
    ''' Lagrange polynomial for extrapolation 
        zeros on xi except for i '''
    res = 1.0
    for s in range(len(xi)):
        if i != s:
            res = res * (x - xi[s]) / (xi[i] - xi[s])
    return res

def ExtrapGen(p):
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
    solPoint = solPointGen(p) # A optimiser : on peut le generer 1! fois AS

    solPoint[0] = -1.0
    solPoint[p+2] = 1.0
    
    ns = p+3
    
    Extrap=np.zeros([ns,ns]);
    for i in range(ns):
        for j in range(ns):
            Extrap[i,j]=lagrange(solPoint[i],solPoint,j); # row = lagrange polynomial non zero on i | column = lagrange pol value on solpoint j
   
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
    somme = 0
    if x in xi:
        ind = np.linspace(1,len(xi),len(xi))
        ind = int(ind[x==xi])
        for s in range(len(xi)):
            if i != s and ind != s:
                res = res * (x - xi[s]) / (xi[i] - xi[s])
                #somme = somme + 1/(x -xi[s])
        #der = somme*res
        der = res / (x[s]-x[ind]) # AS -> il me semble qu'il y avait une erreur avec la "somme" ou alors des computations inutiles !
        
    else:
        for s in range(len(xi)):
            if i != s:
                res = res * (x - xi[s]) / (xi[i] - xi[s])
                somme = somme + 1/(x -xi[s])
        der = somme*res
        
    return der

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

''' Part 4 Derivation of the flux'''

#def DGen(p, phi=1.0):
#    ''' Compute d/dx discretization with '''
#    ''' the Spectral Difference method order p '''
#    ''' phi = 1.0 <=> upwind  flux '''
#    ''' phi = 0.0 <=> centred flux '''
#
#    ns = p+2
#    solPoint = solPointGen(p)
#    solPoint[0] = -1.0
#    solPoint[p+2] = 1.0
##RP
##    nf = p+2
##    fluxPoint = fluxPointGen(p)
#
#    ''' Compute the derivatives matrix '''
#    #RP D = np.zeros([ns, nf])
#    D = np.zeros([ns, ns]) #RP ns + 1 ?
#    for i in range(ns):
#        for ii in range(nf):
#            #RP D[i, ii] = lagrangeDerivative(solPoint[i], ii, fluxPoint)
#            D[i, ii] = lagrangeDerivative(solPoint[i], ii, solPoint) #RP
#
#
#    return D

def D2Gen(p, phi=1.0):    #RP
    ''' Compute d/dx discretization with '''
    ''' the Spectral Difference method order p '''
    ''' phi = 1.0 <=> upwind  flux '''
    ''' phi = 0.0 <=> centred flux '''

    ns = p+3
    solPoint = solPointGen(p)
    solPoint[0] = -1.0
    solPoint[p+2] = 1.0

    ''' Compute the derivatives matrix '''
    D = np.zeros([ns, ns])
    for i in range(ns):
        for j in range(ns):
            D[i, j] = lagrangeDerivative(solPoint[i], j, solPoint) # row = der lagrange polynomial non zero on i | column = der lagrange pol value on solpoint j


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
