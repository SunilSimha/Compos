# This file calculate the 2-point correlation function#

from __future__ import division
import numpy as np
import matterps


def twopcf(x):
    def func(k):
        f = matterps.normalizedmp(k) / (2 * np.pi ** 2) * k * np.sin(k * x) / x
        # f = matterps.normalizedmp(k) * np.cos(k * x)
        return f
    # xi = integ.quad(func,0,float('+inf'),epsabs = 0.0000001,limit = 500)
    xi = np.zeros(np.size(x))
    
    k = np.linspace(0.0001, 1000, 100000)
    p = matterps.normalizedmp(k)
    k0 = k[0]
    for i in range(1, np.size(k)):
        dk = k[i]-k0
        xi = xi + p[i] / (2 * np.pi ** 2) * k[i] * np.sin(k[i] * x) / x * dk
        k0 = k[i]
    return xi
