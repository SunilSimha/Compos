#This script gives examples of growth factor with respect to a with different set of (w_0, w_1)

import numpy as np
import matplotlib.pyplot as mp
import time
import sys

from compos import const, growthfactor

t0 = time.time()

const.initializecosmo(.3, 0, 0.67, 2.7255/2.7, omq = 0.7, sigma8 = 0.8288)

a = np.linspace(0.0001,1,500)

#growthfactor.generateDE(-1,0)
g1 = np.zeros(np.size(a))
for i in range(0,np.size(g1)):
    g1[i] = growthfactor.growfunc_a(a[i])

const.initializecosmo(0.3, 0, 0.67, 2.7255/2.7, omq = 0.7, sigma8 = 0.8288 , w_0 = -0.8, w_1 = 0,n_s = 0.96)

#growthfactor.generateDE(-0.8,0)
g2 = np.zeros(np.size(a))
for i in range(0,np.size(g1)):
    g2[i] = growthfactor.growfunc_a(a[i])

const.initializecosmo(0.3,0, 0.67, 2.7255/2.7, omq = 0.7, sigma8 = 0.8288 , w_0 = -0.8, w_1 = -0.6 ,n_s = 0.96)

#growthfactor.generateDE(-0.8,-0.6)
g3 = np.zeros(np.size(a))
for i in range(0,np.size(g1)):
    g3[i] = growthfactor.growfunc_a(a[i])

const.initializecosmo(0.3, 0, 0.67, 2.7255/2.7, omq = 0.7, sigma8 = 0.8288 , w_0 = -0.8, w_1 = 0.6,n_s = 0.96)

#growthfactor.generateDE(-0.8,0.6)
g4 = np.zeros(np.size(a))
for i in range(0,np.size(g1)):
    g4[i] = growthfactor.growfunc_a(a[i])

t1 = time.time()

print 'The program lasts' 
print t1-t0
print 'seconds.'

mp.subplot(122)
mp.plot(a[3:],g1[3:]/a[3:]/g1[-1],label = '$w_0 = -1, w_1 = 0$')
mp.plot(a[3:],g2[3:]/a[3:]/g2[-1],label = '$w_0 = -0.8, w_1 = 0$')
mp.plot(a[3:],g3[3:]/a[3:]/g3[-1],label = '$w_0 = -0.8, w_1 = -0.6$')
mp.plot(a[3:],g4[3:]/a[3:]/g4[-1],label = '$w_0 = -0.8, w_1 = 0.6$')
mp.xlabel('a')
mp.ylabel('Growth$(\delta/a)/(\delta/a)_{z = 0}$',fontsize = 'small')
mp.legend(fontsize='small')

mp.subplot(121)
mp.plot(a[3:],g1[3:]/a[3:]/(g1[3]/a[3]),label = '$w_0 = -1, w_1 = 0$')
mp.plot(a[3:],g2[3:]/a[3:]/(g2[3]/a[3]),label = '$w_0 = -0.8, w_1 = 0$')
mp.plot(a[3:],g3[3:]/a[3:]/(g3[3]/a[3]),label = '$w_0 = -0.8, w_1 = -0.6$')
mp.plot(a[3:],g4[3:]/a[3:]/(g4[3]/a[3]),label = '$w_0 = -0.8, w_1 = 0.6$')
mp.xlabel('a')
mp.ylabel('Growth$(\delta/a)/(\delta/a)_{a = 0}$',fontsize = 'small')
mp.legend(fontsize='small',loc=3)
mp.savefig('../../results/growthfactor/growthfactor.jpg')
mp.savefig('../../results/growthfactor/growthfactor.pdf')
mp.show()
