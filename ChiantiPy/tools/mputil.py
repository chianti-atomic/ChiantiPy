'''
functions needed for standard Python multiprocessing module mspectrum
'''
import chianti

def doFfQ(inQ, outQ):
    ''' 
    multiprocessing helper for freefree
    '''
    for inputs in iter(inQ.get, 'STOP'):
        ionS = inputs[0]
        temperature = inputs[1]
        wavelength = inputs[2]
        abund = inputs[3]
        em = inputs[4]
        ff = chianti.core.continuum(ionS, temperature, abundance=abund, em=em)
        ff.freeFree(wavelength)
        outQ.put(ff.FreeFree)
    return
    #
    # ----------------------------------------------
    #
def doFbQ(inQ, outQ):
    ''' 
    multiprocessing helper for freeBound
    '''
    for inputs in iter(inQ.get, 'STOP'):
        ionS = inputs[0]
        temperature = inputs[1]
        wavelength = inputs[2]
        abund = inputs[3]
        em = inputs[4]
        fb = chianti.core.continuum(ionS, temperature, abundance=abund, em=em)
        fb.freeBound(wavelength)
        outQ.put(fb.FreeBound)
    return
    #
    # ----------------------------------------------
    #
def doIonQ(inQueue, outQueue):
    ''' 
    multiprocessing helper for ion, also does two-photon
    '''
    for inpts in iter(inQueue.get, 'STOP'):
        ionS = inpts[0]
        temperature = inpts[1]
        density = inpts[2]
        wavelength = inpts[3]
        wvlRange = [wavelength.min(), wavelength.max()]
        filter = inpts[4]
        allLines = inpts[5]
        abund = inpts[6]
        em = inpts[7]
        doContinuum = inpts[8]
        thisIon = chianti.core.ion(ionS, temperature, density, abundance=abund)
        thisIon.intensity(wvlRange = wvlRange, allLines = allLines, em=em)
        if 'errorMessage' not in sorted(thisIon.Intensity.keys()):
            thisIon.spectrum(wavelength,  filter=filter)
        outList = [ionS, thisIon]
        if not thisIon.Dielectronic and doContinuum:
            if (thisIon.Z - thisIon.Ion) in [0, 1]:
                thisIon.twoPhoton(wavelength)
                outList.append(thisIon.TwoPhoton)
        outQueue.put(outList)
    return
