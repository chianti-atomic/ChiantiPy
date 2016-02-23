'''
functions needed for IPython multiprocessing module ipymspectrum
'''
from IPython.parallel import require
import copy
import chianti
@require('copy', 'chianti')
def doFf(inpt):
    ''' multiprocessing helper for freefree'''
    ionS = inpt[0]
    temperature = inpt[1]
    wavelength = inpt[2]
    abund = inpt[3]
    cont = chianti.core.continuum(ionS, temperature, abundance=abund)
    cont.freeFree(wavelength)
    return [ionS, copy.deepcopy(cont.FreeFree)]
    #
    # ----------------------------------------------
    #
@require('copy', 'chianti')
def doFb(inpt):
    '''
    multiprocessing helper for freeBound
    '''
    ionS = inpt[0]
    temperature = inpt[1]
    wavelength = inpt[2]
    abund = inpt[3]
    cont = chianti.core.continuum(ionS, temperature, abundance=abund)
    cont.freeBound(wavelength)
    return [ionS, copy.deepcopy(cont)]
    #
    # ----------------------------------------------
    #
@require('copy', 'chianti')
def doIon(inpt):
    '''
    multiprocessing helper for ion, also does two-photon
    '''
 #     [ionS, temperature, eDensity, wavelength, filter, allLines, abund]
    ionS = inpt[0]
    temperature = inpt[1]
    density = inpt[2]
    wavelength = inpt[3]
    wvlRange = [wavelength.min(), wavelength.max()]
    filter = inpt[4]
    allLines = inpt[5]
    abund = inpt[6]
    em = inpt[7]
    doContinuum = inpt[8]
    thisIon = chianti.core.ion(ionS, temperature, density, abundance=abund)
    thisIon.intensity(wvlRange = wvlRange, allLines = allLines, em=em)
    if 'errorMessage' not in thisIon.Intensity.keys():
        thisIon.spectrum(wavelength,  filter=filter, allLines=allLines)
#        outList = [ionS, thisIon.Spectrum]
#    outList = [ionS, thisIon.Spectrum, thisIon.Intensity]
#    outList = [ionS, thisIon]
#    outList = [ionS, ionS]
    outList = [ionS, thisIon]
    if not thisIon.Dielectronic and doContinuum:
        if (thisIon.Z - thisIon.Ion) in [0, 1]:
            thisIon.twoPhoton(wavelength, em=em)
            outList.append(thisIon.TwoPhoton)
    return outList
    #
    # ----------------------------------------------
    #
@require('copy', 'chianti')
def doAll(inpt):
    '''
    to process ff, fb and line inputs
    '''
    ionS = inpt[0]
    calcType = inpt[1]
    if calcType == 'ff':
        temperature = inpt[2]
        wavelength = inpt[3]
        abund = inpt[4]
        cont = chianti.core.continuum(ionS, temperature, abundance=abund)
        cont.freeFree(wavelength)
        return [ionS, calcType, copy.deepcopy(cont.FreeFree)]
    elif calcType == 'fb':
        temperature = inpt[2]
        wavelength = inpt[3]
        abund = inpt[4]
        cont = chianti.core.continuum(ionS, temperature, abundance=abund)
        cont.freeBound(wavelength)
        return [ionS, calcType, copy.deepcopy(cont)]
    elif calcType == 'line':
        temperature = inpt[2]
        density = inpt[3]
        wavelength = inpt[4]
        wvlRange = [wavelength.min(), wavelength.max()]
        filter = inpt[5]
        allLines = inpt[6]
        abund = inpt[7]
        em = inpt[8]
        doContinuum = inpt[9]
        thisIon = chianti.core.ion(ionS, temperature, density, abundance=abund)
        thisIon.intensity(wvlRange = wvlRange, allLines = allLines, em=em)
        if 'errorMessage' not in thisIon.Intensity.keys():
            thisIon.spectrum(wavelength,  filter=filter, allLines=allLines)
    #        outList = [ionS, thisIon.Spectrum]
    #    outList = [ionS, thisIon.Spectrum, thisIon.Intensity]
    #    outList = [ionS, thisIon]
    #    outList = [ionS, ionS]
        outList = [ionS, calcType, copy.deepcopy(thisIon)]
        if not thisIon.Dielectronic and doContinuum:
            if (thisIon.Z - thisIon.Ion) in [0, 1]:
                thisIon.twoPhoton(wavelength, em=em)
                outList.append(thisIon.TwoPhoton)
        return outList
    #
    # ----------------------------------------------
    #
@require('copy', 'chianti')
def doAll1(inpt):
    '''
    to process ff, fb and line inputs
    '''
    ionS = inpt[0]
    calcType = inpt[1]
    if calcType == 'ff':
        temperature = inpt[2]
        wavelength = inpt[3]
        abund = inpt[4]
        cont = chianti.core.continuum(ionS, temperature, abundance=abund)
        cont.freeFree(wavelength)
        return [ionS, calcType, copy.deepcopy(cont.FreeFree)]
    elif calcType == 'fb':
        temperature = inpt[2]
        wavelength = inpt[3]
        abund = inpt[4]
        cont = chianti.core.continuum(ionS, temperature, abundance=abund)
        cont.freeBound(wavelength)
        return [ionS, calcType, copy.deepcopy(cont)]
    elif calcType == 'line':
        temperature = inpt[2]
        density = inpt[3]
        wavelength = inpt[4]
        wvlRange = [wavelength.min(), wavelength.max()]
        filter = inpt[5]
        allLines = inpt[6]
        abund = inpt[7]
        em = inpt[8]
        doContinuum = inpt[9]
#        thisIon = chianti.core.ion(ionS, temperature, density, abundance=abund)        
        thisIon = ch.ion(ionS, temperature, density, abundance=abund)
        thisIon.intensity(wvlRange = wvlRange, allLines = allLines, em=em)
        if 'errorMessage' not in thisIon.Intensity.keys():
            thisIon.spectrum(wavelength,  filter=filter, allLines=allLines)
    #        outList = [ionS, thisIon.Spectrum]
    #    outList = [ionS, thisIon.Spectrum, thisIon.Intensity]
    #    outList = [ionS, thisIon]
    #    outList = [ionS, ionS]
        outList = [ionS, calcType, copy.deepcopy(thisIon)]
        if not thisIon.Dielectronic and doContinuum:
            if (thisIon.Z - thisIon.Ion) in [0, 1]:
                thisIon.twoPhoton(wavelength, em=em)
                outList.append(thisIon.TwoPhoton)
        return outList
    
