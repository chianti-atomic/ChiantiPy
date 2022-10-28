"""
Line profile filters for creating synthetic spectra.
"""

import numpy as np


def gaussianR(wvl, wvl0, factor=1000.):
    """
    A gaussian filter where the gaussian width is given by `wvl0`/`factor`.

    Parameters
    -----------
    wvl : `~numpy.ndarray`
        Wavelength array
    wvl0 : `~numpy.float64`
        Wavelength filter should be centered on.
    factor : `~numpy.float64`
        Resolving power
    """
    std = wvl0/factor
    wvl = np.asarray(wvl)
    return np.exp(-0.5*((wvl - wvl0)/std)**2)/(np.sqrt(2.*np.pi)*std)


def gaussian(wvl, wvl0, factor=1.):
    """
    A gaussian filter

    Parameters
    -----------
    wvl : `~numpy.ndarray`
        Wavelength array
    wvl0 : `~numpy.float64`
        Wavelength filter should be centered on.
    factor : `~numpy.float64`
        Gaussian width

    integrated value is unity
    """
    wvl = np.asarray(wvl, np.float64)
    return np.exp(-0.5*((wvl - wvl0)/factor)**2)/(np.sqrt(2.*np.pi)*factor)


def boxcar(wvl, wvl0, factor=None):
    """
    Box-car filter

    Parameters
    -----------
    wvl : `~numpy.ndarray`
        Wavelength array
    wvl0 : `~numpy.float64`
        Wavelength filter should be centered on.
    factor : `~numpy.float64`
        Full width of the box-car filter
    """
    wvl = np.asarray(wvl, np.float64)
    dwvl = wvl - np.roll(wvl, 1)
    dwvl[0] = dwvl[1]
    one = np.ones_like(wvl)
    zed = np.zeros_like(wvl)
    if factor is None:
        factor = dwvl.min()
    if factor < dwvl.min():
       raise ValueError('Width must be at least equal to the wavelength step')
    good1 = (wvl > wvl0 - factor/2.)
    good2 = (wvl < wvl0 + factor/2.)
    realgood = np.logical_and(good1, good2)
    return np.where(realgood, one, zed)/(factor)



def lorentz(wvl, wvl0, factor=1.):
    """
    Lorentz profile filter with the exception that all factors are in wavelength units
    rather than frequency as the lorentz profile is usually defined.

    Parameters
    -----------
    wvl : `~numpy.ndarray`
        Wavelength array
    wvl0 : `~numpy.float64`
        Wavelength filter should be centered on.
    factor : `~numpy.float64`
        Value of the so-called constant gamma

    integrated value is unity
    the FWHM is 2*gamma

    .. math::
        L = \\frac{1}{\pi \gamma} \\frac{ \gamma^2}{(\lambda - \lambda_0)^2 + \gamma^2}

    """
    wvl = np.asarray(wvl)
    sigma = factor
    gamma = sigma*np.sqrt(2.*np.log(2.))
    ltz = (gamma/(np.pi))/((wvl - wvl0)**2 + (gamma)**2)
    return ltz

def moffat(wvl, wvl0, factor=2.5):
    """
    Moffat profile with parameters suited to Chandra Letg spectra

    Parameters
    ----------
    wvl : `~numpy.ndarray`
        Wavelength array
    wvl0 : `~numpy.float64`
        Wavelength the filter is centered on.
    factor : `~numpy.float64`
        Resolving power (TODO: correct description)
    """
    wvl = np.asarray(wvl)
    dwvl = np.abs(wvl[1] - wvl[0])
    moffat = 1./(1.+((wvl - wvl0)/0.0275)**2)**factor
    return moffat/(dwvl*moffat.sum())

def voigt(wvl, wvl0, factor=(0.5, 1.)):
    ''' pseudo-Voigt filter
    the sum of a Gaussian and a Lorentzian

    Parameters
    ----------
    wvl : `~numpy.ndarray`
        Wavelength array
    wvl0 : `~numpy.float64`
        Wavelength the filter is centered on.
    factor: array-type
        contains the following 2 parameters
    A : `~numpy.float64`
        relative size of gaussian and lorentz components
        must be between 0. and 1. but this is not currently checked
    sigma:  `~numpy.float64`
        the gaussian width of the gaussian profile (the standard deviation)
        also creates the lorentz component with the same fwhm
    '''
    A = factor[0]
    sigma = factor[1]
    return A*gaussian(wvl, wvl0, sigma) + (1.-A)*lorentz(wvl, wvl0, sigma)
