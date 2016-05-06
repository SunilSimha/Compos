#This script calls transfunction.py to calculate transfunction with no oscillation.

from __future__ import division
import numpy as np
import matplotlib.pyplot as mp
import sys
import os
import matplotlib.gridspec as gridspec

f = os.getcwd()

from compos import const
import compos.transfunction as tf

const.initializecosmo()

k = np.linspace(-3,2)
k = 10 ** k
t = tf.transfunction(k,const.cosmo)
t_nowiggle = tf.t_nowiggle(k,const.cosmo)

mp.loglog(k / 0.7, t, label = 'Wiggled')
mp.loglog(k / 0.7, t_nowiggle, label = 'Dewiggled')
mp.legend()
plottextl = mp.xlim()[0] * 2
plottextd = mp.ylim()[0] * 2
mp.text(plottextl,plottextd, r'$h = $'+str(const.cosmo['h'])+r'$, \Omega_b h^2 = $'+str(const.cosmo['h'] ** 2 * const.cosmo['omega_b'])+r'$, \Omega_c h^2 =$'+str(const.cosmo['h'] ** 2 * const.cosmo['omega_c']),fontsize = 15)
mp.ylabel('T(k)')
mp.xticks([])
mp.savefig('../../results/transfunction/nowiggle.jpg',dpi = 200) 
#mp.savefig('../../results/transfunction/compare.pdf')
mp.show()

