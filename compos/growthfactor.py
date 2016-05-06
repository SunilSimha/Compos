# This code calculates growth factor at a given z or a. Reference is : 0208087

from __future__ import division
import numpy as np
from scipy.integrate import odeint
import scipy.integrate as integrate
import const

# Return a with respect to z


def z2a(z):
    return 1 / (1 + z)

# Return z with respect to a


def a2z(a):
    return 1 / a - 1

# Return equation of state of DE wrt z


def w_z(z):
    w_0 = const.cosmo['w_0']
    w_1 = const.cosmo['w_1']
    return w_0 + w_1 * z / (1 + z)

# Return equation of state of DE wrt z


def w_a(a):
    z = a2z(a)
    return w_z(z)


def func(x):
    return (1 + w_z(x)) / (1 + x)


def exponent(z):
    return integrate.quad(func, 0, z)[0]

# Return Hubble constant wrt z under DE model


def hubble(z):
    global om0, omk, omq
    om0 = const.cosmo['omega_0']
    omq = const.cosmo['omega_q']
    omk = 1 - om0 - omq
    h = const.cosmo['h']
    H0 = h * 100
    H = H0 * np.sqrt((1 + z) ** 3 * om0 - (1 + z) **
                     2 * omk + np.exp(3 * exponent(z)) * omq)
    return H

# Return the density of matter at redshift z


def om_0(z):
    global omk, omq, om0
    om0 = const.cosmo['omega_0']
    omq = const.cosmo['omega_q']
    omk = 1 - om0 - omq
    h = const.cosmo['h']
    H0 = h * 100
    return om0 * (1 + z) ** 3 * (H0 / hubble(z)) ** 2

# Return the density of DE at redshift z


def om_q(z):
    global omk, omq, om0
    om0 = const.cosmo['omega_0']
    omq = const.cosmo['omega_q']
    omk = 1 - om0 - omq
    h = const.cosmo['h']
    H0 = h * 100
    return omq * np.exp(3 * exponent(z)) * (H0 / hubble(z)) ** 2


def ffunc(y, la):
    global omk, omq, om0
    om0 = const.cosmo['omega_0']
    omq = const.cosmo['omega_q']
    omk = 1 - om0 - omq
    x = a2z(np.exp(la))
    f = y[0]
    return [- f ** 2 - (1 - om_0(x) / 2 -
                        (1 + 3 * w_z(x)) / 2 * om_q(x)) *
            f + 3 / 2 * om_0(x), f]

# growth factor wrt a#


def growfunc_a(a):
    la = np.log(a)
    t = [np.log(0.0001), la]
    f0 = 1/(0.0001-1)+1
    ini = [f0, 1]
    d = np.exp(odeint(ffunc, ini, t)[1][1])
    return d

# growth factor nomalized to z = 0#


def growfunc_z(z):
    return (growfunc_a(z2a(z)) / growfunc_a(1))
