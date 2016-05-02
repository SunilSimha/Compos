#This code defines functions that call 'callcamb.py' to calculate transferfunction, then use the transferfunction to calculate matter power spectrum.

import numpy as np
import os
import matplotlib.pyplot as mp

import cambparam, callcamb, const

#Calculate mass power of perturbation from transferfunction calculated by CAMB at z = 0#

def delta2camb():
    kt = callcamb.readtransfunction()
    h = const.cosmo['h']
    k = kt[0]
    t = kt[6] 
    om0 = const.cosmo['omega_0']
    omb = const.cosmo['omega_b']
    Lambda = 1 - om0
    H0 = h * 100
    oml = Lambda / (1 * H0 **2) #omega_Lambda
    omr = 1 - om0 - oml
    c = 2.998 * 10 ** 8 #the velocity of light
    n = 0.96 #power index, set to 1 for Harrison-Zel'dovich-Peebles scale invariant case
    n_tilde = n - 1
    if (Lambda == 0):
        delta_h = 1.95 * 10 ** (-5) * om0 ** (-0.35 - 0.19 * np.log(om0) - 0.17 * n_tilde) * np.e ** (- n_tilde - 0.14 * n_tilde ** 2)
    else:
        delta_h = 1.94 * 10 ** (-5) * om0 ** (-0.785 - 0.05 * np.log(om0)) * np.e ** (- 0.95 * n_tilde - 0.169 * n_tilde ** 2)
#    delta = 2.1 * 10 ** (-9) * k ** (3 + n) * t ** 2
#    delta_h = (k / 0.05) ** ((n - 1) / 2) * 1.91 * 10 ** (-5)
    delta = delta_h ** 2 * (c * k / H0 / 1000) ** (3 + n) * t ** 2
    kdelta = kt
    kdelta[0] = k
    kdelta[1] = delta
    return kdelta
    
def cambsigma8(): #Read sigma8 from cambscreen file. This sigma8 is calculated by camb directly,
    path = const.cosmo['path']
    h = const.cosmo['h']
    hubble = const.cosmo['h'] * 100
    ombh2 = const.cosmo['omega_b'] * h ** 2
    omch2 = const.cosmo['omega_c'] * h ** 2
    omnuh2 = const.cosmo['omega_nu'] * h ** 2
    omk = const.cosmo['omega_k']
    temp_cmb = const.cosmo['T_CMB'] * 2.7
    cpath =  os.path.join(path,'scripts/results/callcamb/spectra')
    os.chdir(cpath)
    param = open('cambscreen'+str(ombh2)+','+str(omch2)+','+str(omnuh2)+','+str(omk)+','+str(hubble)+','+str(temp_cmb)+'.txt', 'r')
    text = param.read()
    tsigma = text.find('matter)=')
    sigma8 = float(text[tsigma+10:-1])
    return sigma8

def mpowercamb(): #Calculate mass power spectrum given by transferfunction calculated by CAMB
    kdelta = delta2camb()
    csigma8 = cambsigma8()
    def sigma8(): #calculate the normalization of the power spectrum (often use sigma(8))
        k = kdelta[0] 
        sigma2 = 0
        def func(x):
            def j1(z):
                j = (z * np.cos(z) - np.sin(z)) / (z ** 2)
                return j
            return (1 / x) * (3 * j1(8 * x)/(8 * x)) ** 2
        delta = delta2camb()[1]
        for i in range(1, np.size(k)):
            sigma2 = sigma2 + delta[i] * func(k[i]) * (k[i]-k[i-1])
        return sigma2
    kdelta[1] = kdelta[1] * 2 * np.pi ** 2 / (((kdelta[0]) ** 3) ) * csigma8 ** 2 / sigma8()
    return kdelta

def tmpowercamb():
    kt = callcamb.readtransfunction()
    t = kt[6]
    p = 2.1 * 10 ** (-9) * (kt[0] * 0.7 / 0.05) ** 0.96 * t ** 2
    kt[1] = p
    return kt

def readmatterps(): #Read matter power spectrum calculated directly by CAMB.
    path = const.cosmo['path']
    cpath = os.path.join(path,'scripts/results/callcamb/spectra')
    os.chdir(cpath)
    h = const.cosmo['h']
    hubble = const.cosmo['h'] * 100
    ombh2 = const.cosmo['omega_b'] * h ** 2
    omch2 = const.cosmo['omega_c'] * h ** 2
    omnuh2 = const.cosmo['omega_nu'] * h ** 2
    omk = const.cosmo['omega_k']
    temp_cmb = const.cosmo['T_CMB'] * 2.7
    kp = np.loadtxt('test('+str(ombh2)+','+str(omch2)+','+str(omnuh2)+','+str(omk)+','+str(hubble)+','+str(temp_cmb)+')_matterpower.dat')
    kp = np.transpose(kp)
    return kp

def plotmatterps(): #This function plots the matterps function calculated by function mpowercamb().
    kpower = mpowercamb()
    k = kpower[0]
    matterps = kpower[1]
    mp.loglog(k,matterps)
    mp.savefig('powerspectrum('+str(ombh2)+','+str(omch2)+','+str(omnuh2)+','+str(omk)+','+str(hubble)+','+str(temp_cmb)+').pdf')
    mp.show()
    return

def comparematterps(): #This function compares mass power spectrums calculated by my code and CAMB.
    kpower = mpowercamb()
    kpowercamb = readmatterps()
    k1 = kpower[0]
    k2 = kpowercamb[0]
    p1 = kpower[1]
    p2 = kpowercamb[1]
    h = const.cosmo['h']
    hubble = const.cosmo['h'] * 100
    ombh2 = const.cosmo['omega_b'] * h ** 2
    omch2 = const.cosmo['omega_c'] * h ** 2
    omnuh2 = const.cosmo['omega_nu'] * h ** 2
    omk = const.cosmo['omega_k']
    temp_cmb = const.cosmo['T_CMB'] * 2.7
    mp.loglog(k1,p1,label = 'Yan\'s code')
    mp.loglog(k2,p2,label = 'CAMB')
    mp.xlim(0.001,1)
    mp.ylim(1,1000000)
    mp.legend()
    mp.text(0.0015,2, '$h = 0.7, \Omega_bh^2 = 0.0226, \Omega_ch^2 = 0.112$',fontsize = 20)
    mp.xlabel('k (h Mpc$^{-1}$)')
    mp.ylabel('P(k) (h$^{-3}$ Mpc$^{-3}$)')
    path = const.cosmo['path']
    cpath = os.path.join(path,'scripts/results/callcamb/')
    os.chdir(cpath)
    mp.savefig('comparespectrum('+str(ombh2)+','+str(omch2)+','+str(omnuh2)+','+str(omk)+','+str(hubble)+','+str(temp_cmb)+').pdf')
    mp.savefig('comparespectrum('+str(ombh2)+','+str(omch2)+','+str(omnuh2)+','+str(omk)+','+str(hubble)+','+str(temp_cmb)+').jpg')
    mp.show()
    return
    
