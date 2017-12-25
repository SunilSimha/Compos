# Compos
This is a instruction document of COMPOS (COdes for Matter POwer Spectrum) written by Ziang Yan. It is written in Python. This fork (me Sunil Simha) is mainly because the software had gotten out of date and required some cleanup for it to install smoothly.

# Changes
* `setuppackages.py` has been renamed `setup.py` and the previously existing file of that name has been deleted. `setuppath.py` now runs in the renamed file.
* All `import` function calls for packages in  `compos` have been modified to `from . import`
* Minor fixes to accommodate change in `Python` syntax.

This is the version July 2015. (forked on Dec 2017)
Prerequisites:
=============
        
        Python
        Python packages: setuptools, numpy, matplotlib, scipy (use sudo apt-get install python-matplotlib to get matplotlib and numpy; $sudo apt-get install python-scipy to get scipy and $sudo apt-get install python-setuptools to get setuptools package)
        CAMB (if you need to call a camb routine):http://camb.info/

Installing:
===========
        
           Before installation, first edit setpath.py. Type the path of CAMB package (for example: '/home/software/camb')in your computer. If you do not have CAMB, then leave it as "''".
           Install: $python setup.py install

Running:
========
        
        Before calling any function in the COMPOS package, you need to initialize a cosmological model by calling:
        const.initializecosmo(om0, omb, H, T_CMB, omq, omnu = 0, omk = 0, w_0 = -1, w_1 = 0, n_s = 1, sigma8 = 0.8 , z = 0, kmax = 20, nonlinear = 0, callcamb = 0)
        The definition and default value of those parameters are (from Planck 2013):
            om0 :    0.316,     #Total density of mass (scaled with critical density)
            omb :    0.049, 	#Density of baryon
            omq :    0.6825,    #Density of dark energy
            omnu :   0, 	#Density of neutrino
            omk :    0,         #Density of curvature
            H :      0.6711, 	#Hubble constant
            T_CMB :  2.728, 	#Temperature of CMB. Normalized to 2.7K
            w_0 :    -1, 	#DE state eq wrt z as w_0+w_1 * z / (1 + z)             
            w_1 :    0, 
            sigma8 : 0.8344, 	#Total fluctuation amplitude within 8 Mpc h^(-1)
            n_s :    0.96, 	#Power index
            z:       1          #Redshift for initialize a nonlinear calculation routine
       
COMPOS package contains five parts:
==================================           

           transfunction: codes calculating transfer functions at given cosmological parameters.
           powerspec:codes calculating linear matter power spectrum and two point correlation function.
           growthfactor:codes calculating growthfactor at given dark energy (w_0,w_1) model at given a or z.
           callcamb:codes calling CAMB package calculating transfer function and matter power spectrum.
           pyhalofit:python version for HALOFIT (see arxiv:1208.2701).

           Each package has script.py files for examples to use them.
              
