#This code calls 'const.py' to calculate the mass power spectrum.
#Formulas in this code are also from arxiv:9709112v1
#NOTE: all spectrum calculated by this code depends on k scaled with Mpc^-1. If you need to use h Mpc^-1 as unit, please set yourself.

from __future__ import division
import numpy as np
import scipy.integrate as integ
import matplotlib.pyplot as mp

import sys
import os

import const
import transfunction as tf
import growthfactor

#global parameters#

global c,s8
s8 = 0
c = 2.998 * 10 ** 8 #the velocity of light

#density of mass at redshift z#

def omega(z): 
    global h, om0, omb, Lambda, H0, oml, omr, c, s8_given, ns
    h = const.cosmo['h']
    om0 = const.cosmo['omega_0']
    omb = const.cosmo['omega_b']
    s8_given = const.cosmo['sigma8']
    ns = const.cosmo['n_s']
    Lambda = 1 - om0
    H0 = h * 100
    oml = Lambda / (1 * H0 **2) #omega_Lambda
    omr = 1 - om0 - oml
    omega_z = om0 * (1 + z) ** 3 / (oml + omr * (1 + z) ** 2 + om0 * (1 + z) ** 3)
    return omega_z

#density of Lambda at redshift of z#

def omegaL(z): 
    global h, om0, omb, Lambda, H0, oml, omr, c, s8_given, ns
    omega_L = oml / (oml + omr * (1 + z) ** 2 + om0 * (1 + z) ** 3)
    return omega_L

#calculate the power of perturbation

def delta2(k): 
    global h, om0, omb, Lambda, H0, oml, omr, c, s8_given, ns
    h = const.cosmo['h']
    om0 = const.cosmo['omega_0']
    omb = const.cosmo['omega_b']
    s8_given = const.cosmo['sigma8']
    ns = ns = const.cosmo['n_s']
    Lambda = 1 - om0
    H0 = h * 100
    oml = Lambda / (1 * H0 **2) #omega_Lambda
    omr = 1 - om0 - oml
    n = ns #power index, set to 1 for Harrison-Zel'dovich-Peebles scale invariant case
    n_tilde = n - 1
    if (Lambda == 0):
        delta_h = 1.95 * 10 ** (-5) * om0 ** (-0.35 - 0.19 * np.log(om0) - 0.17 * n_tilde) * np.e ** (- n_tilde - 0.14 * n_tilde ** 2)
    else:
        delta_h = 1.94 * 10 ** (-5) * om0 ** (-0.785 - 0.05 * np.log(om0)) * np.e ** (- 0.95 * n_tilde - 0.169 * n_tilde ** 2)
    k_eq = tf.keq(const.cosmo)
    k_silk = tf.ksilk(const.cosmo)
    t = tf.transfunction(k,const.cosmo)

    delta = delta_h ** 2 * (c * k / H0 / 1000) ** (3 + n) * t ** 2
    return delta

#calculate the mass power spectrum

def matterps(k): 
    d = delta2(k)
    p = d * 2 * np.pi ** 2 / (k ** 3)
    return p

#calculate the normalization of the power spectrum (often use sigma(8)). Note that k is in the unit of Mpc^-1#

def sigma8(): 
    global s8
    if s8 == 0:
        def func(x):
            def j1(z):
                j = (z * np.cos(z) - np.sin(z)) / (z ** 2)
                return j
            return (1 / x) * delta2(x) * (3 * j1(8 / const.cosmo['h'] * x)/(8 / const.cosmo['h'] * x)) ** 2
        sigma2 = integ.quad(func, 0, 20)
        s8 = np.sqrt(sigma2[0])
        return
    else:
        return

#calculate normalized power spectrum(by sigma8)#

def normalizedmp(k,z = 0): 
    global s8
    sigma8()
    normalp = matterps(k) * (s8_given / s8) ** 2
    if z == 0:
        return normalp
    else:
        return normalp * (growthfactor.growfunc_z(z) / growthfactor.growfunc_z(0)) ** 2
    
#Delta^2 without wiggle

def delta2nw(k): 
    global h, om0, omb, Lambda, H0, oml, omr, c, s8_given, ns
    h = const.cosmo['h']
    om0 = const.cosmo['omega_0']
    omb = const.cosmo['omega_b']
    s8_given = const.cosmo['sigma8']
    ns = const.cosmo['n_s']
    Lambda = 1 - om0
    H0 = h * 100
    oml = Lambda / (1 * H0 **2) #omega_Lambda
    omr = 1 - om0 - oml
    n = 0.96 #power index, set to 1 for Harrison-Zel'dovich-Peebles scale invariant case
    n_tilde = n - 1
    if (Lambda == 0):
        delta_h = 1.95 * 10 ** (-5) * om0 ** (-0.35 - 0.19 * np.log(om0) - 0.17 * n_tilde) * np.e ** (- n_tilde - 0.14 * n_tilde ** 2)
    else:
        delta_h = 1.94 * 10 ** (-5) * om0 ** (-0.785 - 0.05 * np.log(om0)) * np.e ** (- 0.95 * n_tilde - 0.169 * n_tilde ** 2)
    k_eq = tf.keq(const.cosmo)
    k_silk = tf.ksilk(const.cosmo)
    t = tf.t_nowiggle(k,const.cosmo)
    delta = delta_h ** 2 * (c * k / H0 / 1000) ** (3 + n) * t ** 2
    return delta

#Matter power spectrum without wiggle

def matterpsnw(k): 
    d = delta2nw(k)
    p = d * 2 * np.pi ** 2 / (k ** 3)
    return p

#calculate normalized power spectrum(by sigma8) without wiggle#

def normalizedmpnw(k,z = 0): 
    global s8
    sigma8()
    normalp = matterpsnw(k) * (s8_given / s8) ** 2
    if z == 0:
        return normalp
    else:
        return normalp * (growthfactor.growfunc_z(z) / growthfactor.growfunc_z(0)) ** 2

def kstar():
    sum = integ.quad(normalizedmp,0,3000)[0]
    kstar = ((1 / (3 * np.pi ** 2)) * sum) ** (-0.5)
    return kstar


