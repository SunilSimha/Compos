#This script calculates linear matter power spectrum at z = 0 with a given set of cosmologycal parameters#

import numpy as np
import matplotlib.pyplot as mp

from compos import const, matterps

const.initializecosmo()

k = np.linspace(-4,np.log10(6.71600e+02),10000)
k = 10 ** k * const.cosmo['h']
p = matterps.normalizedmp(k)

#p_nw = matterps.normalizedmpnw(k)

mp.loglog(k / const.cosmo['h'], p * const.cosmo['h'] ** 3)
plottextl = mp.xlim()[0] * 2
plottextd = mp.ylim()[0] * 2
mp.text(plottextl,plottextd, r'$h = $'+str(const.cosmo['h'])+', $\Omega_mh^2$ = '+str(const.cosmo['omega_0'] * const.cosmo['h'] ** 2)+', $\Omega_bh^2$ = '+str(const.cosmo['omega_b'] * const.cosmo['h'] ** 2),fontsize = 15)
mp.ylabel('P(k)(h$^{-3}$Mpc$^{-3}$)')
mp.xlabel('k (h Mpc$^{-1}$)')
mp.savefig('../../results/matterps/mps.jpg') 
mp.savefig('../../results/matterps/mps.pdf')
mp.show()
