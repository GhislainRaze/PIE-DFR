import numpy as np
import matplotlib.pyplot as plt


x0=np.zeros(1000000)
sol0=np.zeros(len(x0))
i=0
file= open('results_init.txt', "r") 
for line in file:
    parts = line.split(",")
    x0[i]=parts[0]
    sol0[i]=parts[1]   
    i=i+1
        
file.close()




x1=np.zeros(1000000)
yth1=np.zeros(len(x1))
sol1=np.zeros(len(x1))
i=0


# Choose order p
p=3
file= open('results_CFL=0.7_p='+str(p)+'.txt', "r") 
for line in file:
    parts = line.split(",")
    
    x1[i]=parts[0]
    yth1[i]=parts[1]
    sol1[i]=parts[2]   
    i=i+1
        
file.close()


x2=np.zeros(len(x1))
yth2=np.zeros(len(x1))
sol2=np.zeros(len(x1))
i=0
file= open('avbp21d_fic.out_3', "r") 
for line in file:
    parts = line.split(" ")
    x2[i]=parts[0]
    sol2[i]=parts[1]   
    i=i+1
        
file.close()

#x4=np.zeros(len(x1))
#yth4=np.zeros(len(x1))
#sol4=np.zeros(len(x1))
#i=0
#file= open('resultsdiff_CFL=0.01_p=2.txt', "r") 
#for line in file:
#    parts = line.split(",")
#    x4[i]=parts[0]
#    yth4[i]=parts[1]
#    sol4[i]=parts[2]   
#    i=i+1
        
#file.close()

#x5=np.zeros(len(x1))
#yth5=np.zeros(len(x1))
#sol5=np.zeros(len(x1))
#i=0
#file= open('resultsdiff_CFL=0.1_p=6.txt', "r") 
#for line in file:
#    parts = line.split(",")
#    x5[i]=parts[0]
#    yth5[i]=parts[1]
#    sol5[i]=parts[2]   
#    i=i+1
        
#file.close()

#x6=np.zeros(len(x1))
#yth6=np.zeros(len(x1))
#sol6=np.zeros(len(x1))
#i=0
#file= open('results_CFL=0.001_p=6.txt', "r") 
#for line in file:
#    parts = line.split(",")
#    x6[i]=parts[0]
#    yth6[i]=parts[1]
#    sol6[i]=parts[2]   
#    i=i+1
        
#file.close()
plt.plot(x0,sol0, 'k-')
plt.plot(x1,yth1, 'b--')
plt.plot(x1,sol1, 'c-')
plt.plot(x2,sol2, 'r-')
#plt.plot(x4,sol4, 'r-')
#plt.plot(x6,sol6, 'g-')


plt.legend(['initial','exact','SD p='+str(p),'AVBP' ])
plt.xlabel('x(m)')
plt.ylabel('Yfic')

plt.show()


