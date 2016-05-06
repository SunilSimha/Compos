# This script uses twopcf.py to calculate a 2-point
# correlation function for a given set of cosmological parameters#

import numpy as np
import matplotlib.pyplot as mp
import time
from compos import matterps, twopcf, const

t0 = time.time()

const.initializecosmo()

r = np.linspace(1, 300, 1000)
xi = twopcf.twopcf(r)

k = np.linspace(-4, 1, 10000)
k = 10 ** k
p = matterps.normalizedmp(k)

mp.subplot(121)
mp.loglog(k / 0.7, p * 0.7 ** 3)
plottextl = mp.xlim()[0] * 2
plottextd = mp.ylim()[0] * 2
mp.text(plottextl, plottextd,
        '$h = 0.7, \Omega_bh^2 = 0.0226, \Omega_ch^2 = 0.112$', fontsize=15)
mp.ylabel('$P(k)(h^{-3}$ $Mpc^{-3}$)')
mp.xlabel('$k(h$ $Mpc^{-1}$)')

mp.subplot(122)
mp.plot(r * 0.7, xi * (r * 0.7) ** 2)
mp.xlabel('$r(Mpc/h)$')
mp.ylabel(r'$r^2\xi (Mpc/h)$')
t1 = time.time()
print 'The program lasts ' + str(t1 - t0) + ' seconds'
mp.savefig('../../results/matterps/twopcf.jpg', dpi=200)
mp.show()


