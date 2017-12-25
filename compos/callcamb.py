# This routine calls camb to calculate transferfunction with
# a given set of cosmological parameters.

import numpy as np
import os
from . import const
from . import cambparam

# write into a param.ini file.


global h, hubble, ombh2, omch2, omnuh2, omk, temp_cmb, path, cambpath
path = const.cosmo['path']
cambpath = const.cosmo['cambpath']


def initialize(kmax=10):
    cambparam.writeparam(const.cosmo, kmax)
    cambparam.writeoutput(const.cosmo)
    return

# call a CAMB routine#


def cambbypy(kmax=10):
    global path, cambpath
    initialize(kmax)
    h = const.cosmo['h']
    hubble = const.cosmo['h'] * 100
    ombh2 = const.cosmo['omega_b'] * h ** 2
    omch2 = const.cosmo['omega_c'] * h ** 2
    omnuh2 = const.cosmo['omega_nu'] * h ** 2
    omk = const.cosmo['omega_k']
    temp_cmb = const.cosmo['T_CMB'] * 2.7
    os.chdir(cambpath)
    os.path.join(path, 'scripts/results/callcamb/spectra')
    os.system('./camb paramsformps.ini > ' + path +
              '/scripts/results/callcamb/spectra/cambscreen' +
              str(ombh2) + ',' + str(omch2) + ',' + str(omnuh2) +
              ',' + str(omk) + ',' + str(hubble) + ',' + str(temp_cmb) +
              '.txt')
    return

# read transferfunction from .dat file.#


def readtransfunction():
    global path
    path = const.cosmo['path']
    cambbypy()
    h = const.cosmo['h']
    hubble = const.cosmo['h'] * 100
    ombh2 = const.cosmo['omega_b'] * h ** 2
    omch2 = const.cosmo['omega_c'] * h ** 2
    omnuh2 = const.cosmo['omega_nu'] * h ** 2
    omk = const.cosmo['omega_k']
    temp_cmb = const.cosmo['T_CMB'] * 2.7
    os.path.join(path, 'scripts/results/callcamb/spectra')
    os.chdir(path)
    kt = np.loadtxt('test(' + str(ombh2) + ',' + str(omch2) +
                    ',' + str(omnuh2) + ',' + str(omk) + ',' +
                    str(hubble) + ',' + str(temp_cmb) + ')_transfer_out.dat')
    kt = np.transpose(kt)
    return kt

# Read matter power spectrum calculated directly by CAMB.


def readmatterps():
    path = const.cosmo['path']
    cpath = os.path.join(path, 'scripts/results/callcamb/spectra')
    os.chdir(cpath)
    h = const.cosmo['h']
    hubble = const.cosmo['h'] * 100
    ombh2 = const.cosmo['omega_b'] * h ** 2
    omch2 = const.cosmo['omega_c'] * h ** 2
    omnuh2 = const.cosmo['omega_nu'] * h ** 2
    omk = const.cosmo['omega_k']
    temp_cmb = const.cosmo['T_CMB'] * 2.7
    kp = np.loadtxt('test(' + str(ombh2) + ',' + str(omch2) + ',' +
                    str(omnuh2) + ',' + str(omk) + ',' + str(hubble) +
                    ','+str(temp_cmb)+')_matterpower.dat')
    kp = np.transpose(kp)
    return kp

# Read sigma_8 calculated directly by CAMB.


def readsigma8():
    global path
    path = const.cosmo['path']
    cambbypy()
    h = const.cosmo['h']
    hubble = const.cosmo['h'] * 100
    ombh2 = const.cosmo['omega_b'] * h ** 2
    omch2 = const.cosmo['omega_c'] * h ** 2
    omnuh2 = const.cosmo['omega_nu'] * h ** 2
    omk = const.cosmo['omega_k']
    temp_cmb = const.cosmo['T_CMB'] * 2.7
    os.path.join(path, 'scripts/results/callcamb/spectra')
    os.chdir(path)
    param = open('cambscreen' + str(ombh2) + ',' + str(omch2) + ',' +
                 str(omnuh2) + ',' + str(omk) + ',' + str(hubble) + ',' +
                 str(temp_cmb) + '.txt', 'r')
    text = param.read()
    sig = text.find('r)=')
    sigma8 = float(text[sig + 3:-1])
    return sigma8
    
