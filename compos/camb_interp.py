# This file read a matter power spectrum calculated by a CAMB
# routine and do interpolation to get P_k. Then de-wiggle it.

import numpy as np
from scipy import interpolate
import scipy.integrate as integ

from . import callcamb
from . import growthfactor
from . import const

f = interpolate.UnivariateSpline([1, 2], [1, 2], k=1)

# Notice: the k and P(k) in a CAMB output file is in the unit of
# h Mpc^-1 and h^-3 Mpc ^3.Here I first unify the units with matterps.py#

global s8
s8 = 0


def inimp(kmax=10):
    global f, k, s8
    callcamb.cambbypy(kmax)
    kp = callcamb.readmatterps()
    k = kp
    k[0] = kp[0] * const.cosmo['h']
    k[1] = kp[1] / (const.cosmo['h'] ** 3)
    f = interpolate.UnivariateSpline(k[0], k[1], k=1, s=0)
    return

# give an unnormalized interpolated power spectrum function
# based on the camb output .dat file. Note that the input k is with unit Mpc^-1


def CAMB_Pk(k):
    global f, s8
    p = f(k)
    return p


def sigma8():
    global s8

    def func(x):
        def j1(z):
            j = (z * np.cos(z) - np.sin(z)) / (z ** 2)
            return j
        return (1 / x) * CAMB_delta2(x) * (3 * j1(8 / const.cosmo['h'] * x) /
                                           (8 / const.cosmo['h'] * x)) ** 2
    sigma2 = integ.quad(func, 0, 10)
    s8 = np.sqrt(sigma2[0])
    return s8

# give a normalized interpolated power spectrum function
# based on the camb output .dat file. Note that the input k is with unit Mpc^-1


def normCAMB_Pk(k, z=0):
    global f, s8
    s8 = sigma8()
    p = CAMB_Pk(k) * (const.cosmo['sigma8'] / s8) ** 2
    if z == 0:
        return p
    else:
        return p * (growthfactor.growfunc_z(z) /
                    growthfactor.growfunc_z(0)) ** 2


def CAMB_delta2(k):
    return CAMB_Pk(k) * k ** 3 / (2 * np.pi ** 2)


def kstar():
    global s8
    s8 = sigma8()
#    callcamb.cambbypy(kmax = 2000)
    kp = callcamb.readmatterps()
    k = kp
    k[0] = kp[0]  # * const.cosmo['h']
    k[1] = kp[1] * (const.cosmo['sigma8'] / s8) ** 2
    # / (const.cosmo['h'] ** 3)
    k0 = k[0][0]
    sum = 0
    for i in range(1, np.size(k[0])-1):
        dk = k[0][i] - k0
        sum = sum + k[1][i] * dk
        k0 = k[0][i]
    kstar = ((1 / (3 * np.pi ** 2)) * sum) ** (-0.5)
    delta2 = k[1] * k[0] ** 3 / (2 * np.pi ** 2)
    
    def j1(z):
        j = (z * np.cos(z) - np.sin(z)) / (z ** 2)
        return j
    sigma2 = 0
    for i in range(1, np.size(k[0])-1):
        sigma2 = sigma2 + (1 / k[0][i]) * \
                 delta2[i] * (3 * j1((8) * k[0][i]) /
                              ((8) * k[0][i])) ** 2 * (k[0][i] - k[0][i-1])
    s_8 = np.sqrt(sigma2)
    print(s_8)
    return kstar
