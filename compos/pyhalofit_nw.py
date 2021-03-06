# This code is a python version of halofit.
# It gives nonlinear correction on dewiggled matter power spectrum.
# The formulas are in Appendix of 1208.2701.
# The dewiggle routine is given by 9709112.

from __future__ import division
import numpy as np
from . import matterps
from . import const
from . import growthfactor


def Sigma2(R, z=0, gf=1):
    if matterps.s8 == 0:
        matterps.sigma8()
    nint = 3000
    sum1 = 0.0
    sum2 = 0.0
    sum3 = 0.0
    if (z == 0 and gf == 1):
        d = 1
    elif (z != 0 and gf == 1):
        d = growthfactor.growfunc_z(z) / growthfactor.growfunc_z(0)
    elif (z == 0 and gf != 1):
        d = gf
    for i in range(1, nint):
        t = (float(i) - 0.5) / float(nint)
        k = -1 + 1.0 / t
        y = k
        dt2 = matterps.normalizedmpnw(k) * d ** 2 * k ** 3 / (2 * np.pi ** 2)
        x = k * R
        w1 = np.exp(- x * x)
        w2 = 2 * x * x * w1
        w3 = 4 * x * x * (1 - x * x) * w1
        sum1 = sum1 + w1 * dt2 / y / t / t
        sum2 = sum2 + w2 * dt2 / y / t / t
        sum3 = sum3 + w3 * dt2 / y / t / t
    sum1 = sum1 / float(nint)
    sum2 = sum2 / float(nint)
    sum3 = sum3 / float(nint)
    d1 = -sum2 / sum1
    d2 = -sum2 * sum2 / sum1 / sum1 - sum3 / sum1
    return np.array([sum1, d1, d2])
        

def solveforr(rd=0, d=1):  # Solve for k_sigma^(-1) with sigma^2(R) = 1
    global om0, Rnl, neff, C, an, bn, cn
    global gamman, alphan, betan, mun, nun, f1, f2, f3
    i = 1
    xlogr1 = -5.0
    xlogr2 = 5
    while i == 1:
        rmid = (xlogr2 + xlogr1) / 2.0
        rmid = 10 ** rmid
        r = Sigma2(rmid, z=rd, gf=d)
        diff = r[0] - 1.0
        if (abs(diff) <= 0.001):
            Rnl = rmid
            neff = -3-r[1]
            C = -r[2]
            return
        elif (diff > 0.001):
            xlogr1 = np.log10(rmid)
            continue
        else:
            xlogr2 = np.log10(rmid)
            continue


def supparams(z=0, growthf=1):  # supplementary parameters(A5-A13)
    global om0, Rnl, neff, C, an, bn, cn, gamman
    global alphan, betan, mun, nun, f1, f2, f3
    solveforr(rd=z, d=growthf)
    w = growthfactor.w_z(z)
    om0 = const.cosmo['omega_0']
    if (w == -1):
        omq = 0
    else:
        omq = growthfactor.omq(z)
    an = 10 ** (1.5222 + 2.8553 * neff + 2.3706 * neff **
                2 + 0.9903 * neff ** 3 + 0.225 * neff **
                4 - 0.6038 * C + 0.1749 * omq * (1 + w))
    bn = 10 ** (-0.5642 + 0.5864 * neff + 0.5716 * neff **
                2 - 1.5474 * C + 0.2279 * omq * (1 + w))
    cn = 10 ** (0.3698 + 2.0404 * neff + 0.8161 * neff ** 2 + 0.5869 * C)
    gamman = 0.1971 - 0.0843 * neff + 0.8460 * C
    alphan = np.abs(6.0835 + 1.3373 * neff - 0.1959 * neff ** 2 - 5.5274 * C)
    betan = 2.0379 - 0.7354 * neff + 0.3157 * neff ** 2 + \
        1.249 * neff ** 3 + 0.3980 * neff ** 4 - 0.1682 * C
    mun = 0
    nun = 10 ** (5.2105 + 3.6902 * neff)
    om_q = const.cosmo['omega_q']
    if(abs(1-om0) > 0.01):
        f1a = om0 ** (-0.0732)
        f2a = om0 ** (-0.1423)
        f3a = om0 ** (0.0725)
        f1b = om0 ** (-0.0307)
        f2b = om0 ** (-0.0585)
        f3b = om0 ** (0.0743)
        frac = om_q/(1.0-om0)
        f1 = frac * f1b + (1 - frac) * f1a
        f2 = frac * f2b + (1 - frac) * f2a
        f3 = frac * f3b + (1 - frac) * f3a
    else:
        f1 = 1.0
        f2 = 1.0
        f3 = 1.0
    return


def deltaq(k, z=0):  # Two halo terms(A2)
    global om0, Rnl, neff, C, an, bn, cn
    global gamman, alphan, betan, mun, nun, f1, f2, f3
    y = k * Rnl

    def funy(x):
        return x / 4 + x ** 2 / 8
    deltal = matterps.normalizedmpnw(k, z) * k ** 3 / (2 * np.pi ** 2)
    deltaq = deltal * ((1 + deltal) ** betan / (1 + alphan * deltal)) * \
        np.exp(- funy(y))
    return deltaq


def deltah(k):  # One halo terms(A3)
    global om0, Rnl, neff, C, an, bn, cn
    global gamman, alphan, betan, mun, nun, f1, f2, f3
    y = k * Rnl
    deltaprime = an * y ** (3 * f1) / (1 + bn * y **
                                       f2 + (cn * f3 * y) ** (3 - gamman))
    deltah = deltaprime / (1 + mun * y ** (-1) + nun * y ** (-2))
    return deltah


def nldelta_nw(k, z=0):  # Total perturbation
    return deltaq(k, z) + deltah(k)


def nlpowerspec_nw(k, z=0):  # Spectrum of nonlinear spectrum
    return nldelta_nw(k, z) * k ** (-3) * 2 * np.pi ** 2
