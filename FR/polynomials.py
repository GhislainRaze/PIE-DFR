# -*- coding: utf-8 -*-
"""
Created on Sun Oct 30 16:10:19 2016

@author: Ghislain Raze
"""

import numpy as np
import matplotlib.pyplot as plt


#############################################################################################################
#############################################################################################################
##                                                                                                         ##
##                                                                                                         ##
##                                          POLYNOMIALS.PY                                                 ##
##                       Polynomials used for the Flux Reconstruction approach                             ##
##                                                                                                         ##
##                                                                                                         ##
#############################################################################################################
#############################################################################################################




# Returns the coefficients of the Legendre polynomials of given order
# The coefficients are arranged in a matrix such that
#       The rows contain the Legendre polynomial coefficients of the same order
#       The columns correspond to the coefficient of a given power of x
def legendreCoef(order):
    
    if order < 0 or order-int(order)!=0 :
        return 0.
    elif order == 0:
        return 1.
    elif order == 1:
        return np.array([[1,0],[0,1]])
    else:
        # Construction from recurrence relation        
        leg = np.zeros([order+1,order+1])
        leg[0][0] = 1
        leg[1][1] = 1
        for k in range(2,order+1):
            for l in range(k+1):
                leg[k][l] = (2*k-1.)/k * leg[k-1][l-1] - (k-1.)/k * leg[k-2][l]
        return leg  


# Returns the Legendre polynomials evaluated at x
# The returned matrix is such that the lines corresponds to the Legendre polynomials of the same order evaluated at x
# startOrder can be used to only compute the order-startOrder+1 last polynomials
def legendre(x,order,startOrder = 0):

    y = np.zeros([order+1-startOrder,len(x)])
    leg = legendreCoef(order)
    for i in range(order+1):
        for j in range(startOrder,order+1):
            y[j-startOrder] = y[j-startOrder] + leg[j,i]*x**i

    return y


# Returns the derivatives of the Legendre polynomials evaluated at x
# The returned matrix is such that the lines corresponds to the derivatives of the Legendre polynomials of the same order evaluated at x
# startOrder can be used to only compute the order-startOrder+1 last polynomials
def legendreDerivatives(x,order,startOrder = 0):

    y = np.zeros([order+1-startOrder,len(x)])
    leg = legendreCoef(order)
    for i in range(1,order+1):
        for j in range(startOrder,order+1):
            y[j-startOrder] = y[j-startOrder] + i*leg[j,i]*x**(i-1)

    return y

# Returns the coefficients of the Radau polynomials of given order 
# The coefficients are arranged in a matrix such that
#       The rows contain the Radau polynomial coefficients of the same order
#       The columns correspond to the coefficient of a given power of x       
def radauCoef(order):
    leg = legendreCoef(order)
    rad = np.zeros([order+1,order+1])
    for k in range(1,order+1):
        rad[k] = 0.5*(-1.)**k * (leg[k] - leg[k-1])
        
    return rad
 

# Returns the Radau polynomials evaluated at x
# The returned matrix is such that the lines corresponds to the Radau polynomials of the same order evaluated at x
# startOrder can be used to only compute the order-startOrder+1 last polynomials
def radau(x,order,startOrder = 0):

    y = np.zeros([order+1-startOrder,len(x)])
    rad = radauCoef(order)
    for i in range(order+1):
        for j in range(startOrder,order+1):
            y[j-startOrder] = y[j-startOrder] + rad[j,i]*x**i

    return y


# Returns the derivatives of the Radau polynomials evaluated at x
# The returned matrix is such that the lines corresponds to the derivatives of the Radau polynomials of the same order evaluated at x
# startOrder can be used to only compute the order-startOrder+1 last polynomials
def radauDerivatives(x,order,startOrder = 0):

    y = np.zeros([order+1-startOrder,len(x)])
    rad = radauCoef(order)
    for i in range(1,order+1):
        for j in range(startOrder,order+1):
            y[j-startOrder] = y[j-startOrder] + i*rad[j,i]*x**(i-1)
    return y



# Returns the coefficients of the Lobatto polynomials of given order
# The coefficients are arranged in a matrix such that
#       The rows contain the Lobatto polynomial coefficients of the same order
#       The columns correspond to the coefficient of a given power of x    
def lobattoCoef(order):
    leg = legendreCoef(order)
    lob = np.zeros([order+1,order+1])
    for k in range(2,order+1):
        lob[k] = (leg[k] - leg[k-2])
        
    return lob


# Returns the Lobatto polynomials evaluated at x
# The returned matrix is such that the lines corresponds to the lobatto polynomials of the same order evaluated at x
# startOrder can be used to only compute the order-startOrder+1 last polynomials  
def lobatto(x,order,startOrder = 0):

    y = np.zeros([order+1-startOrder,len(x)])
    lob = lobattoCoef(order)
    for i in range(order+1):
        for j in range(startOrder,order+1):
            y[j-startOrder] = y[j-startOrder] + lob[j,i]*x**i

    return y    
    
    
# Returns the derivatives of the Lobatto polynomials evaluated at x
# The returned matrix is such that the lines corresponds to the derivatives of the lobatto polynomials of the same order evaluated at x
# startOrder can be used to only compute the order-startOrder+1 last polynomials  
def lobattoDerivatives(x,order,startOrder = 0):

    y = np.zeros([order+1-startOrder,len(x)])
    lob = lobattoCoef(order)
    for i in range(1,order+1):
        for j in range(startOrder,order+1):
            y[j-startOrder] = y[j-startOrder] + i*lob[j,i]*x**(i-1)

    return y       


# Returns the values of the Lagrange polynomials of given order (useless function)   
#def lagrange(x,xi,yi):
#    
#    lag = np.ones([len(xi),len(x)])     
#    dx = np.ones([len(xi),len(x)])
#    
#    for i in range(len(dx)):
#        dx[i] = x - xi[i]
#        
#    for i in range(len(lag)):
#        for j in range(len(lag)):
#            if i != j:
#                lag[i] = lag[i]*dx[j]/(xi[i]-xi[j])
#        lag[i] = lag[i]*yi[i]
#    
#    y = np.sum(lag,1)
#    
#    return y
    

# Returns the g2 function evaluated at x
def g2(x,p):
    rad = radau(x,p,p-1)
    y = (p-1)/(2*p-1.0) *rad[1] + p/(2*p-1.0)*rad[0]
    return y


# Returns the derivatives of g2 function evaluated at x
def g2Derivative(x,p):
    drad = radauDerivatives(x,p,p-1)
    y = (p-1)/(2*p-1.0) *drad[1] + p/(2*p-1.0)*drad[0]
    return y
    


# Returns the gl and gr functions evaluated at x
def glgr(x,p):
    p = p-1
    leg = legendre(x,p,p)
    gl = 0.5*(-1)**p * (1-x) * leg[0]
    gr = 0.5 * (1+x) * leg[0]
    return gl, gr
  
# Returns the derivatives of the gl and gr functions evaluated at x
def glgrDerivatives(x,p):
    p = p-1
    leg = legendre(x,p,p)
    dleg = legendreDerivatives(x,p,p)
    gl = 0.5*(-1)**p * ( (1-x) *dleg[0] - leg[0] )
    gr = 0.5 * ( (1+x) * dleg[0] + leg[0] )
    return gl, gr    


# Energy-preserving schemes ?






# To check the correction functions
"""
p = 3
x = np.linspace(-1,1,100)
dx = x[1]-x[0]

gl, gr = glgr(x,p)


plt.figure()
plt.plot(x,gl)
plt.plot(x,gr)

dglNum = np.zeros(len(x))

for i in range(1,len(x)-1):
    dglNum[i] = 0.5*(gl[i+1]-gl[i-1])/dx

dglNum[0] = (gl[1]-gl[0])/dx
dglNum[-1] = (gl[-1]-gl[-2])/dx

dgrNum = np.zeros(len(x))

for i in range(1,len(x)-1):
    dgrNum[i] = 0.5*(gr[i+1]-gr[i-1])/dx

dgrNum[0] = (gr[1]-gr[0])/dx
dgrNum[-1] = (gr[-1]-gr[-2])/dx

dgl,dgr = glgrDerivatives(x,p)
plt.figure()
plt.plot(x,dgl)
plt.plot(x,dglNum,'r--')


plt.figure()
plt.plot(x,dgr)
plt.plot(x,dgrNum,'r--')


gl = g2(x,p)
gr = np.flipud(gl)
dgl = g2Derivative(x,p)
dgr = -np.flipud(dgl)


plt.figure()
plt.plot(x,gl)
plt.plot(x,gr)

dglNum = np.zeros(len(x))

for i in range(1,len(x)-1):
    dglNum[i] = 0.5*(gl[i+1]-gl[i-1])/dx

dglNum[0] = (gl[1]-gl[0])/dx
dglNum[-1] = (gl[-1]-gl[-2])/dx

dgrNum = np.zeros(len(x))

for i in range(1,len(x)-1):
    dgrNum[i] = 0.5*(gr[i+1]-gr[i-1])/dx

dgrNum[0] = (gr[1]-gr[0])/dx
dgrNum[-1] = (gr[-1]-gr[-2])/dx

plt.figure()
plt.plot(x,dgl)
plt.plot(x,dglNum,'r--')


plt.figure()
plt.plot(x,dgr)
plt.plot(x,dgrNum,'r--')
plt.show()
"""