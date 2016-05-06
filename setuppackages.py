# This is the setup file of compos codes#
# Run it with 'sudo python setup.py install'#
# The code should be re-installed after any change#

from setuptools import setup, find_packages

packages = find_packages()
setup(
    name="Compos",
    version="0",
    packages=packages,
    install_requires=['numpy', 'scipy', 'matplotlib', ],

    #    author = "Ziang Yan",
    #    author_email = "yanzastro@163.com",
    #    description = "Codes for Matter Power Spectrum",
    #    keywords = ("astronomy cosmology cosmological matter power spectrum"),
    #    license = "THCA",
    #    long_description = \
    # """
    # Compos contains five parts:
    #
    #    transfunction: codes calculating transfer functions at
    #    given cosmological parameters.
    #    powerspec:codes calculating linear matter power spectrum
    #    and two point correlation function.
    #    growthfactor:codes calculating growthfactor at given dark
    #    energy (w_0,w_1) model at given a or z.
    #    callcamb:codes calling CAMB package calculating transfer
    #    function and matter power spectrum.
    #    pyhalofit:python version for halofit.
    #
    # Each package has script.py files for examples to use them.
    # """
)


