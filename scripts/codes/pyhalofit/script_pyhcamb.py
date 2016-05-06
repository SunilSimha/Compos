import numpy as np
import matplotlib.pyplot as mp

import os
from compos import const, camb_interp, pyhalofit_camb

f = os.getcwd()
const.initializecosmo(callcamb=1)

pyhalofit_camb.supparams()
x = np.linspace(-4, np.log(2), 1000)
k = 10 ** x
p = pyhalofit_camb.nlpowerspec(k)
t = np.zeros((2, np.size(k)))
t[0] = k
t[1] = p
t = np.transpose(t)
p2 = camb_interp.CAMB_Pk(k)

os.chdir(f)
np.savetxt('../../results/pyhalofit/halofit_camb.txt', t)
mp.loglog(k, p, label='nonlinear')
mp.loglog(k, p2, label='linear')
mp.text(0.001, 0.1, '$h = 0.7, \Omega_mh^2$ = ' +
        str(const.cosmo['omega_0']*const.cosmo['h']**2) +
        ',$ \Omega_qh^2 = $'+str(const.cosmo['omega_q'] *
                                 const.cosmo['h']**2)+',', fontsize=15)
mp.legend()
mp.xlabel('k (h Mpc$^{-1}$)')
mp.ylabel('P(k) (h$^{-3}$ Mpc$^{-3}$)')
mp.savefig('../../results/pyhalofit/pyhalofit_camb.pdf')
mp.savefig('../../results/pyhalofit/pyhalofit_camb.jpg')
mp.show()
