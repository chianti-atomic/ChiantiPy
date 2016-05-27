''' line profile filters from creating synthetic spectra

Copyright 2009, 2010 Kenneth P. Dere

This software is distributed under the terms of the GNU General Public License
that is found in the LICENSE file

'''
import numpy as np
def gaussianR(wvl,wvl0, factor=1000.):
    '''
    a gaussian filter where factor is the resolving power, so that the gaussian width (standard deviation)
    is given by wvl0/factor'''
    if factor:
        std = wvl0/factor
    else:
        print(' the resolving power of the gaussianR filter is undefined')
        return None
    wvl = np.asarray(wvl, 'float64')
    return np.exp(-0.5*((wvl - wvl0)/std)**2)/(np.sqrt(2.*np.pi)*std)

def gaussian(wvl,wvl0, factor=0):
    '''
    a gaussian filter where factor is the gaussian width (standard deviation)
    '''
    if factor:
        std = factor
    else:
        print(' the width of the gaussian filter is undefined')
        return None
    wvl = np.asarray(wvl, 'float64')
    dwvl = wvl - np.roll(wvl, 1)
    dwvl[0] = dwvl[1]
    return np.exp(-0.5*((wvl - wvl0)/std)**2)/(np.sqrt(2.*np.pi)*std)
    #
def boxcar(wvl, wvl0, factor=0):
    ''' box-car filter, factor is the full width of the box filter'''
    wvl = np.asarray(wvl, 'float64')
    dwvl = wvl - np.roll(wvl, 1)
    dwvl[0] = dwvl[1]
    one = np.ones_like(wvl)
    zed = np.zeros_like(wvl)
    if factor:
        # width must be at least equal to the wavelength step
        width = max(factor, dwvl.min())
        print((' width = %10.2e'%(width)))
    else:
        print(' the width of the box filter is undefined')
        return None
    good1 = (wvl > wvl0 - width/2.)
    good2 = (wvl < wvl0 + width/2.)
    realgood = np.logical_and(good1, good2)
    return np.where(realgood, one, zed)/(width)
    #
def lorentz(wvl, wvl0, factor=0):
    '''the lorentz profile with the exception that all factors are in wavelength units
    rather than frequency as the lorentz profile is usually defined.
    factor is the value of the so-called constant gamma'''
    if factor:
        gamma = factor
    else:
        print(' the factor gamma of the lorentz filter is undefined')
        return None
    wvl = np.asarray(wvl, 'float64')
    dwvl = wvl - np.roll(wvl, 1)
    dwvl[0] = dwvl[1]
    ltz = (gamma/(2.*np.pi)**2)/((wvl - wvl0)**2 + (gamma/(4.*np.pi))**2)
    return np.abs(ltz/(dwvl*ltz.sum()))
    #
def moffat(wvl, wvl0, factor=2.5):
    '''
    the moffat profile with parameters suited to Chandra Letg spectra
    '''
    wvl = np.asarray(wvl, 'float64')
    dwvl = np.abs(wvl[1] - wvl[0])
    moffat = 1./(1.+((wvl - wvl0)/0.0275)**2)**factor
    return moffat/(dwvl*moffat.sum())
