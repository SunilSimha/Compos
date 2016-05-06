#This file calls transfunction.py to calculate transfunction for a set of given cosmological parameters.

from __future__ import division
import numpy as np
import matplotlib.pyplot as mp
import sys
import os

from compos import const
import compos.transfunction as tf

const.initializecosmo()

k = np.linspace(-3,2)
k = 10 ** k
k = k * const.cosmo['h']
t = tf.transfunction(k,const.cosmo)

ob = const.cosmo['omega_b']
o0 = const.cosmo['omega_0']


mp.loglog(k/const.cosmo['h'],np.abs(t),label = 'COMPOS')
#mp.ylim(0.001,10)

pl = mp.xlim()[0] * 2
pd = mp.ylim()[0] * 2
pu = mp.ylim()[0] * 6

mp.xlabel(r'$k(h Mpc^{-1}$)')
mp.ylabel(r'$|T(k)|$')
mp.text(pl,pu,'$\Omega_b=$'+str(ob))
mp.text(pl,pd,'$\Omega_0=$'+str(o0)+', h='+str(const.cosmo['h']))
mp.legend(fontsize = 10)

mp.savefig('../../results/transfunction/transf(h=' + str(const.cosmo['h']) + ',O_0=' + str(const.cosmo['omega_0']) + ',O_b=' + str(const.cosmo['omega_b']) + ').pdf')
mp.show()

