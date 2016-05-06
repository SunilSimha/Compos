# This code contains main constant or paramefters
# for cosmology and an initialize function.

# Cosmological parameters#

global cosmo
cosmo = {'omega_0': 0.316, 	# Total density of mass
         'omega_b': 0.049, 	# Density of baryon
         'omega_c': 0.267, 	# Density of CDM
         'omega_q': 0.6825, 	# Density of DE
         'omega_nu': 0, 	# Density of neutrino
         'omega_k': 0, 	# Density of curvature
         'h': 0.11, 		# Hubble constant
         'T_CMB': 2.728/2.7, 	# Temperature of CMB. Normalized to 2.7K
         'w_0': -1, 		# DE state eq wrt z as w_0+w_1 * z / (1 + z)
         'w_1': 0,
         'sigma8': 0.8344, 	# Total fluctuation amplitude within 8 Mpc h^(-1)
         'n_s': 1, 		# Power index
         'b_g': 1.1, 		# Bias on MPS
         'path': '/home/yanza15/research/software/Compos',  # path of compos
         'cambpath': '/home/yanza15/research/software/camb/'  # path of camb
         }

# Initialize a cosmological model


def initializecosmo(om0=0.316, omb=0.049, H=67.11, T_CMB=2.728, omq=0.6825,
                    omnu=0, w_0=-1, w_1=0, n_s=0.96, sigma8=0.8344, b_g=1,
                    z=0, kmax=20, nonlinear=0, callcamb=0):
    # kmax is the maxmine of k in a camb MPS.
    # For nonlinear and callcamb, 1 means generate
    # nonlinear model or initialize a camb routine.
    global cosmo
    cosmo['omega_0'] = om0
    cosmo['omega_b'] = omb
    cosmo['omega_c'] = om0 - omb
    cosmo['omega_q'] = omq
    cosmo['omega_nu'] = omnu
    cosmo['omega_k'] = 1 - omq - om0 - omnu
    cosmo['h'] = H / 100
    cosmo['T_CMB'] = T_CMB / 2.7
    cosmo['w_0'] = w_0
    cosmo['w_1'] = w_1
    cosmo['n_s'] = n_s
    cosmo['sigma8'] = sigma8
    cosmo['b_g'] = b_g
    import pyhalofit
    import camb_interp
    if nonlinear == 1:
        # Initialize a nonlinear matter power spectrum calculation routine
        pyhalofit.supparams(z)
    if callcamb == 1:
        camb_interp.inimp(kmax=kmax)  # Initialize a camb routine
    return
