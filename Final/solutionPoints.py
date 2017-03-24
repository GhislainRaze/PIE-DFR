# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 07:37:41 2016

@author: Ghislain
"""

import numpy as np

def writeSP(p,SPchoice,z1=0.):

    if(SPchoice == 1):          # Gauss-Lobatto points

        solPoint = np.zeros(p+1)
        for i in range(len(solPoint)):
            solPoint[i] = - np.cos(np.pi * (2. * (i + 1) - 1) / (2 * (p + 1)))  

    elif(SPchoice==2):          # Gauss points

        order = p+1
        if order==1 :
            solPoint=np.array([0.0])
        elif order==2:
            solPoint=np.array([-np.sqrt(1./3),np.sqrt(1./3)])
        elif order==3:
            solPoint=np.array([-0.77459667, 0.0, 0.77459667])
        elif order==4:
            solPoint=np.array([-0.86113631, -0.33998104, 0.33998104, 0.86113631])
        elif order==5:
            solPoint=np.array([-0.90617985, -0.53846931, 0.0, 0.53846931, 0.90617985])
        elif order==6:
            solPoint=np.array([-0.93246951, -0.66120939, -0.23861918, 0.23861918, 0.66120939, 0.93246951])
        elif order==7:
            solPoint=np.array([-0.94910791, -0.74153119, -0.40584515, 0.0, 0.40584515, 0.74153119, 0.94910791])
        elif order==8:
            solPoint=np.array([-0.96028986, -0.79666648, -0.52553241, -0.18343464, 0.18343464, 0.52553241, 0.79666648, 0.96028986])
        else:
            solPoint=0
            print 'Choose different value for number of Gauss integration points'

    elif(SPchoice==3):
        if p==1 :
            solPoint= np.array([-np.sqrt(1./3),np.sqrt(1./3)])
        elif p == 2:
            solPoint= np.array([-np.sqrt(3./5),0.,np.sqrt(3./5)])
        elif p == 3:
            if z1 < 0. or z1 > 1./np.sqrt(3.):
                solPoint= np.array([-np.sqrt(3./7.+2./7.*np.sqrt(6./5.)),-np.sqrt(3./7.-2./7.*np.sqrt(6./5.)),np.sqrt(3./7.-2./7.*np.sqrt(6./5.)),np.sqrt(3./7.+2./7.*np.sqrt(6./5.))])
            else:
                solPoint= np.array([-np.sqrt((3-5*z1**2)/(5-15*z1**2)),-z1,z1,np.sqrt((3-5*z1**2)/(5-15*z1**2))])

        elif p == 4: 
            if z1 < 0. or z1 > np.sqrt(3./5.):
                solPoint= np.array([-1./3.*np.sqrt(5.+2.*np.sqrt(10./7.)),-1./3.*np.sqrt(5.-2.*np.sqrt(10./7.)),0.,1./3.*np.sqrt(5.-2.*np.sqrt(10./7.)),1./3.*np.sqrt(5.+2.*np.sqrt(10./7.))])
            else:
                solPoint= np.array([-np.sqrt(3./7.)*np.sqrt((7*z1**2-5.)/(5*z1**2-3.)),-z1,0.,z1,np.sqrt(3./7.)*np.sqrt((7*z1**2-5.)/(5*z1**2-3.))])
        elif p == 5:
            if z1 < 0. or z1 > np.sqrt(3./7.-2./7.*np.sqrt(6./5.)):
                solPoint = np.array([-0.93246951,-0.66120939,-0.23861919,0.23861919,0.66120939,0.093246951])
            else:
                a = (63.-630.*z1**2+735.*z1**4)
                b = (-70.+588.*z1**2-630.*z1**4)
                c = (15.-70.*z1**2+63.*z1**4)
                delta = b**2-4*a*c
                roots = np.sort(np.array([-np.sqrt((-b-np.sqrt(delta))/(2*a)),np.sqrt((-b-np.sqrt(delta))/(2*a)),-np.sqrt((-b+np.sqrt(delta))/(2*a)),np.sqrt((-b+np.sqrt(delta))/(2*a))]))
                fz1 = roots[2]
                solPoint = np.array([-np.sqrt(1./7.*(6.-7.*z1**2+(1.-42.*z1**2-49.*z1**4)/(-6.+7.*z1**2+63.*z1**4+fz1**2))),-fz1,-z1,z1,fz1,np.sqrt(1./7.*(6.-7.*z1**2+(1.-42.*z1**2-49.*z1**4)/(-6.+7.*z1**2+63.*z1**4+fz1**2)))])

    # Write initialization
    np.savetxt("SP.txt",solPoint)   