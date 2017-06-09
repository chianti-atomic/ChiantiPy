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


def gaussian(wvl, wvl0, factor=1):
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
    """
    wvl = np.asarray(wvl, 'float64')
    dwvl = wvl - np.roll(wvl, 1)
    dwvl[0] = dwvl[1]
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
    wvl = np.asarray(wvl, 'float64')
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


def lorentz(wvl, wvl0, factor=1):
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
    """
    wvl = np.asarray(wvl)
    dwvl = wvl - np.roll(wvl, 1)
    dwvl[0] = dwvl[1]
    ltz = (factor/(2.*np.pi)**2)/((wvl - wvl0)**2 + (factor/(4.*np.pi))**2)
    return np.abs(ltz/(dwvl*ltz.sum()))


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
