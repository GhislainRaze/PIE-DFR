# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 07:37:41 2016

@author: Ghislain
"""

import numpy as np

    
def gaussRule(order):
    
    if order==1 :
        t=np.array([0.0]);
        w=np.array([2.0]);
    elif order==2:
        t=np.array([-np.sqrt(1./3),np.sqrt(1./3)])
        w=np.array([1,1])
    elif order==3:
        t=np.array([-0.77459667, 0.0, 0.77459667])
        w=np.array([0.55555555,0.88888889,0.55555555])
    elif order==4:
        t=np.array([-0.86113631, -0.33998104, 0.33998104, 0.86113631])
        w=np.array([0.34785485,0.65214515,0.65214515,0.34785485])
    elif order==5:
        t=np.array([-0.90617985, -0.53846931, 0.0, 0.53846931, 0.90617985])
        w=np.array([0.23692689,0.47862867,0.56888889,0.47862867,0.23692689])
    elif order==6:
        t=np.array([-0.93246951, -0.66120939, -0.23861918, 0.23861918, 0.66120939, 0.93246951])
        w=np.array([0.17132449,0.36076157,0.46791393,0.46791393,0.36076157,0.17132449])
    elif order==7:
        t=np.array([-0.94910791, -0.74153119, -0.40584515, 0.0, 0.40584515, 0.74153119, 0.94910791])
        w=np.array([0.12948497,0.27970539,0.38183005,0.41795918,0.38183005,0.27970539,0.12948497 ])
    elif order==8:
        t=np.array([-0.96028986, -0.79666648, -0.52553241, -0.18343464, 0.18343464, 0.52553241, 0.79666648, 0.96028986])
        w=np.array([0.10122854,0.22238103,0.31370665,0.36268378,0.36268378,0.31370665,0.22238103,0.10122854])
    else:
        t=0
        w=0
        print 'Choose different value for number of Gauss integration points'
        
    return t,w


def gaussIntegration(sol2,intLagrange,wG,dx):
    integral = 0.
    for i in range(len(sol2)):
        integral = integral + 0.5*dx[i]*np.dot(np.dot(intLagrange,sol2[i,:]),wG)

    return integral 