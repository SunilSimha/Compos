# This code calls camb to generate a matter power spectrum
# at a given set of cosmological parameters and interpolate it at a given k#

# This script calculates linear matter power spectrum at
# z = 0 with a given set of cosmologycal parameters#

import numpy as np
import os
import matplotlib.pyplot as mp

from compos import const, camb_interp

f = os.getcwd()
const.initializecosmo(kmax=20, callcamb=1)
os.chdir(f)

k = np.linspace(-3, 1, 10000)
k = 10 ** k * const.cosmo['h']

p = camb_interp.normCAMB_Pk(k)

mp.loglog(k / const.cosmo['h'], p * const.cosmo['h'] ** 3)
plottextl = mp.xlim()[0] * 2
plottextd = mp.ylim()[0] * 2
mp.text(plottextl, plottextd, r'$h = $' + str(const.cosmo['h']) +
        ', $\Omega_mh^2$ = ' + str(const.cosmo['omega_0'] *
                                   const.cosmo['h'] ** 2) +
        ', $\Omega_bh^2$ = ' + str(const.cosmo['omega_b'] *
                                   const.cosmo['h'] ** 2), fontsize=15)
mp.ylabel('P(k)(h$^{-3}$Mpc$^{-3}$)')
mp.xlabel('k (h Mpc$^{-1}$)')
mp.savefig('../../results/callcamb/cambmps.jpg')
mp.savefig('../../results/callcamb/cambmps.pdf')
mp.show()
