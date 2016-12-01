from datetime import datetime
import copy
import multiprocessing as mp

import numpy as np

import ChiantiPy.tools.data as chdata
import ChiantiPy.tools.constants as const
import ChiantiPy.tools.filters as chfilters
import ChiantiPy.tools.util as util
import ChiantiPy.Gui as chgui
from ChiantiPy.base import ionTrails
from ChiantiPy.base import specTrails
import ChiantiPy.tools.mputil as mputil

defaults = chdata.Defaults


class mspectrum(ionTrails, specTrails):
    ''' this is the multiprocessing version of spectrum
    set proc to the desired number of processors, default=3

    Calculate the emission spectrum as a function of temperature and density.

    temperature and density can be arrays but, unless the size of either is one (1),
    the two must have the same size

    the returned spectrum will be convolved with a filter of the specified width on the
    specified wavelength array

    the default filter is gaussianR with a resolving power of 100.  Other filters,
    such as gaussian, box and lorentz, are available in chianti.filters.  When using the box filter,
    the width should equal the wavelength interval to keep the units of the continuum and line
    spectrum the same.

    A selection of elements can be make with elementList a list containing the names of elements
    that are desired to be included, e.g., ['fe','ni']

    A selection of ions can be make with ionList containing the names of
    the desired lines in Chianti notation, i.e. C VI = c_6

    Both elementList and ionList can not be specified at the same time

    a minimum abundance can be specified so that the calculation can be speeded up by excluding
    elements with a low abundance. With solar photospheric abundances -

    setting minAbund = 1.e-4 will include H, He, C, O, Ne
    setting minAbund = 2.e-5 adds  N, Mg, Si, S, Fe
    setting minAbund = 1.e-6 adds  Na, Al, Ar, Ca, Ni

    Setting em will multiply the spectrum at each temperature by the value of em.

    em [for emission measure], can be a float or an array of the same length as the
    temperature/density.
    allLines = 1 will include lines with either theoretical or observed wavelengths.  allLines=0 will
    include only those lines with observed wavelengths

    proc = the number of processors to use
    timeout - a small but non-zero value seems to be necessary
    '''
    def __init__(self, temperature, eDensity, wavelength, filter=(chfilters.gaussianR, 1000.), label=0, elementList = 0, ionList = 0, minAbund=1.e-6, keepIons=0, abundanceName=0,  doContinuum=1, allLines = 1, em =0,  proc=3, verbose = 0,  timeout=0.1):
        #
        t1 = datetime.now()
        # creates Intensity dict from first ion calculated
        setupIntensity = 0
        #
        #masterlist = chdata.MasterList
        # use the ionList but make sure the ions are in the database
        #if elementList:
            #for i,  one in enumerate(elementList):
                #elementList[i] = one.lower()
            #alist = []
            #for one in masterlist:
                #stuff = util.convertName(one)
                #if stuff['Element'] in  elementList:
                    #alist.append(one)
            #masterlist = alist
        #elif ionList:
            #alist=[]
            #for one in ionList:
                #if masterlist.count(one):
                    #alist.append(one)
                #else:
                    #if verbose:
                        #pstring = ' %s not in CHIANTI database'%(one)
                        #print(pstring)
            #masterlist = alist
        self.Defaults = defaults
        #
#        masterlist = chdata.MasterList
        self.Defaults = defaults
        self.Temperature = np.asarray(temperature, 'float64')
        nTemp = self.Temperature.size
        self.EDensity = np.asarray(eDensity, 'float64')
        nDen = self.EDensity.size
        nTempDen = max([nTemp, nDen])
        self.NTempDen = nTempDen
        #
        if type(em) == int and em == 0:
            em = np.ones(self.NTempDen, 'float64')
        elif type(em) == float and em > 0.:
            em = np.ones(self.NTempDen, 'float64')*em
        elif type(em) == list or type(em) == tuple:
            em = np.asarray(em, 'float64')
        self.Em = em
        #
        #
        if self.Em.any() > 0.:
            ylabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1} \AA^{-1}$ $'
        else:
            ylabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1} \AA^{-1}$ ($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'
        #
        xlabel = 'Wavelength ('+self.Defaults['wavelength'] +')'
        #
        #
        self.AllLines = allLines
        #
        if not abundanceName:
            self.AbundanceName = self.Defaults['abundfile']
        else:
            if abundanceName in chdata.Abundance:
                self.AbundanceName = abundanceName
            else:
                abundChoices = list(chdata.Abundance.keys())
#                for one in wvl[topLines]:
#                    wvlChoices.append('%12.3f'%(one))
                abundChoice = chgui.gui.selectorDialog(abundChoices,label='Select Abundance name')
                abundChoice_idx = abundChoice.selectedIndex
                self.AbundanceName = abundChoices[abundChoice_idx[0]]
                abundanceName = self.AbundanceName
                print(' Abundance chosen:  %s '%(self.AbundanceName))
        #
        abundAll = chdata.Abundance[self.AbundanceName]['abundance']
        self.AbundAll = abundAll
        #
#        nonzed = abundAll > 0.
#        minAbundAll = abundAll[nonzed].min()
#        # if minAbund is even set
#        if minAbund:
#            if minAbund < minAbundAll:
#                minAbund = minAbundAll
        #ionInfo = chio.masterListInfo()
        wavelength = np.asarray(wavelength)
        nWvl = wavelength.size
        self.Wavelength = wavelength
#        wvlRange = [wavelength.min(), wavelength.max()]
        #
        proc = min([proc, mp.cpu_count()])
        #
        freeFree = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        freeBound = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        twoPhoton = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        lineSpectrum = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        #
        #  free-free multiprocessing setup
        ffWorkerQ = mp.Queue()
        ffDoneQ = mp.Queue()
        #
        #  free-bound multiprocessing setup
        #
        fbWorkerQ = mp.Queue()
        fbDoneQ = mp.Queue()
        #
        #  ion multiprocessing setup
        ionWorkerQ = mp.Queue()
        ionDoneQ = mp.Queue()
        #
        self.IonsCalculated = []
        if keepIons:
            self.IonInstances = {}
        self.Finished = []
        #

#        self.Todo = []
        self.ionGate(elementList = elementList, ionList = ionList, minAbund=minAbund, doContinuum=doContinuum, verbose = verbose)
        #
        for akey in sorted(self.Todo.keys()):
            zStuff = util.convertName(akey)
            Z = zStuff['Z']
            abundance = chdata.Abundance[self.AbundanceName]['abundance'][Z - 1]
            if verbose:
                print(' %5i %5s abundance = %10.2e '%(Z, const.El[Z-1],  abundance))
            if verbose:
                print(' doing ion %s for the following processes %s'%(akey, self.Todo[akey]))
            if 'ff' in self.Todo[akey]:
#                if verbose:
#                    print(' doing ff')
                ffWorkerQ.put((akey, temperature, wavelength, abundance, em))
#                allInpt.append([akey, 'ff', temperature, wavelength, abundance])
            if 'fb' in self.Todo[akey]:
#                if verbose:
#                    print(' doing fb')
                fbWorkerQ.put((akey, temperature, wavelength, abundance, em))
#                allInpt.append([akey, 'fb', temperature, wavelength, abundance])
            if 'line' in self.Todo[akey]:
#                if verbose:
#                    print(' doing line')
                ionWorkerQ.put((akey, temperature, eDensity, wavelength, filter, allLines, abundance, em, doContinuum))
#                allInpt.append([akey, 'line', temperature, eDensity, wavelength, filter, allLines, abundance, em, doContinuum])


        #
        ffWorkerQSize = ffWorkerQ.qsize()
        fbWorkerQSize = fbWorkerQ.qsize()
        ionWorkerQSize = ionWorkerQ.qsize()
        if doContinuum:
            ffProcesses = []
            for i in range(proc):
                #-kpd
                p = mp.Process(target=mputil.doFfQ, args=(ffWorkerQ, ffDoneQ))
#                p = mp.Process(target=doFfQ, args=(ffWorkerQ, ffDoneQ))
                p.start()
                ffProcesses.append(p)
    #       timeout is not necessary
            for p in ffProcesses:
                if p.is_alive():
                    p.join(timeout=timeout)
#            for i in range(proc):
#                ffProcesses.append('STOP')
            #
            for iff in range(ffWorkerQSize):
                thisFreeFree = ffDoneQ.get()
                if 'rate' in sorted(thisFreeFree.keys()):
                    freeFree += thisFreeFree['rate']
            for p in ffProcesses:
                if not isinstance(p, str):
                    p.terminate()
        #
            fbProcesses = []
            for i in range(proc):
                #-kpd
                p = mp.Process(target=mputil.doFbQ, args=(fbWorkerQ, fbDoneQ))
#                p = mp.Process(target=doFbQ, args=(fbWorkerQ, fbDoneQ))
                p.start()
                fbProcesses.append(p)
    #       timeout is not necessary
            for p in fbProcesses:
                if p.is_alive():
                    p.join(timeout=timeout)
#            for i in range(proc):
#                fbProcesses.append('STOP')
            #
            for ifb in range(fbWorkerQSize):
                thisFreeBound = fbDoneQ.get()
                if 'rate' in sorted(thisFreeBound.keys()):
                    freeBound += thisFreeBound['rate']
            for p in fbProcesses:
                if not isinstance(p, str):
                    p.terminate()
        #
        ionProcesses = []
        if ionWorkerQSize < proc:
            proc = ionWorkerQSize
        for i in range(proc):
            #-kpd
            p = mp.Process(target=mputil.doIonQ, args=(ionWorkerQ, ionDoneQ))
#            p = mp.Process(target=doIonQ, args=(ionWorkerQ, ionDoneQ))
            p.start()
            ionProcesses.append(p)
#            ionWorkerQ.put('STOP')
#       timeout is not necessary
        for p in ionProcesses:
#            print' process is alive:  ', p.is_alive()
            if p.is_alive():
#                p.join()
                p.join(timeout=timeout)
#        for i in range(proc):
#            ionProcesses.append('STOP')
        #
        for ijk in range(ionWorkerQSize):
            out = ionDoneQ.get()
            ions = out[0]
            if verbose:
                print(' collecting calculation for %s'%(ions))
            thisIon = out[1]
#            thisSpectrum = thisIon.Spectrum
            thisIntensity = thisIon.Intensity
            if not 'errorMessage' in sorted(thisIntensity.keys()):
                self.Finished.append(ions)
                if keepIons:
                    self.IonInstances[ions] = copy.deepcopy(thisIon)
                if setupIntensity:
                    for akey in sorted(self.Intensity.keys()):
                        self.Intensity[akey] = np.hstack((copy.copy(self.Intensity[akey]), thisIntensity[akey]))
                else:
                    setupIntensity = 1
                    self.Intensity  = thisIntensity
                #
                if not 'errorMessage' in sorted(thisIon.Spectrum.keys()):
                    lineSpectrum += thisIon.Spectrum['intensity']
#                if nTempDen == 1:
#                    lineSpectrum += thisSpectrum['intensity']
#                else:
#                    for iTempDen in range(nTempDen):
#                        lineSpectrum[iTempDen] += thisSpectrum['intensity'][iTempDen]
               # check for two-photon emission
                if len(out) == 3:
                    tp = out[2]
                    twoPhoton += tp['rate']
#                    if nTempDen == 1:
#                        twoPhoton += tp['rate']*em[0]
#                    else:
#                        for iTempDen in range(nTempDen):
#                            twoPhoton[iTempDen] += tp['rate'][iTempDen]*em[iTempDen]
            else:
                if 'errorMessage' in sorted(thisIntensity.keys()):
                    print(thisIntensity['errorMessage'])
        #
        for p in ionProcesses:
            if not isinstance(p, str):
                p.terminate()
        #
        #
        #
        self.FreeFree = {'wavelength':wavelength, 'intensity':freeFree.squeeze()}
        self.FreeBound = {'wavelength':wavelength, 'intensity':freeBound.squeeze()}
        self.LineSpectrum = {'wavelength':wavelength, 'intensity':lineSpectrum.squeeze()}
        self.TwoPhoton = {'wavelength':wavelength, 'intensity':twoPhoton.squeeze()}
        #
        total = freeFree + freeBound + lineSpectrum + twoPhoton
        #
        t2 = datetime.now()
        dt=t2-t1
        print(' elapsed seconds = %12.3f'%(dt.seconds))
        #
        if nTempDen == 1:
            integrated = total
        else:
            integrated = total.sum(axis=0)
        #
        if type(label) == type(''):
            if hasattr(self, 'Spectrum'):
                self.Spectrum[label] = {'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'ions':self.IonsCalculated, 'em':em,
                'ions':self.IonsCalculated, 'Abundance':self.AbundanceName, 'xlabel':xlabel, 'ylabel':ylabel}
            else:
                self.Spectrum = {label:{'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'ions':self.IonsCalculated, 'em':em,
                'ions':self.IonsCalculated, 'Abundance':self.AbundanceName, 'xlabel':xlabel, 'ylabel':ylabel}}
        else:
            self.Spectrum ={'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'ions':self.IonsCalculated,
            'em':em, 'Abundance':self.AbundanceName, 'xlabel':xlabel, 'ylabel':ylabel}
    #
    # -------------------------------------------------------------------------
    #
    #-kpd
#def doFfQ(inQ, outQ):
#    """
#    Multiprocessing helper for `ChiantiPy.core.continuum.freeFree`
#
#    Parameters
#    -----------
#    inQ : `~multiprocessing.Queue`
#        Ion free-free emission jobs queued up by multiprocessing module
#    outQ : `~multiprocessing.Queue`
#        Finished free-free emission jobs
#    """
#    for inputs in iter(inQ.get, 'STOP'):
#        ionS = inputs[0]
#        temperature = inputs[1]
#        wavelength = inputs[2]
#        abund = inputs[3]
#        em = inputs[4]
#        #-kpd
##        ff = ch.continuum(ionS, temperature, abundance=abund, em=em)
#        ff = ChiantiPy.core.continuum(ionS, temperature, abundance=abund, em=em)
#        ff.freeFree(wavelength)
#        outQ.put(ff.FreeFree)
#    return
#
#
#def doFbQ(inQ, outQ):
#    """
#    Multiprocessing helper for `ChiantiPy.core.continuum.freeBound`
#
#    Parameters
#    -----------
#    inQ : `~multiprocessing.Queue`
#        Ion free-bound emission jobs queued up by multiprocessing module
#    outQ : `~multiprocessing.Queue`
#        Finished free-bound emission jobs
#    """
#    for inputs in iter(inQ.get, 'STOP'):
#        ionS = inputs[0]
#        temperature = inputs[1]
#        wavelength = inputs[2]
#        abund = inputs[3]
#        em = inputs[4]
#        #-kpd
##        fb = ch.continuum(ionS, temperature, abundance=abund, em=em)
#        fb = ChiantiPy.core.continuum(ionS, temperature, abundance=abund, em=em)
#        fb.freeBound(wavelength)
#        outQ.put(fb.FreeBound)
#    return
#
#
#def doIonQ(inQueue, outQueue):
#    """
#    Multiprocessing helper for `ChiantiPy.core.ion` and `ChiantiPy.core.ion.twoPhoton`
#
#    Parameters
#    -----------
#    inQueue : `~multiprocessing.Queue`
#        Jobs queued up by multiprocessing module
#    outQueue : `~multiprocessing.Queue`
#        Finished jobs
#    """
#    for inpts in iter(inQueue.get, 'STOP'):
#        ionS = inpts[0]
#        temperature = inpts[1]
#        density = inpts[2]
#        wavelength = inpts[3]
#        wvlRange = [wavelength.min(), wavelength.max()]
#        filter = inpts[4]
#        allLines = inpts[5]
#        abund = inpts[6]
#        em = inpts[7]
#        doContinuum = inpts[8]
#        #-kpd
##        thisIon = ch.ion(ionS, temperature, density, abundance=abund)
#        thisIon = ChiantiPy.core.ion(ionS, temperature, density, abundance=abund)
#        thisIon.intensity(wvlRange = wvlRange, allLines = allLines, em=em)
#        if 'errorMessage' not in sorted(thisIon.Intensity.keys()):
#            thisIon.spectrum(wavelength,  filter=filter)
#        outList = [ionS, thisIon]
#        if not thisIon.Dielectronic and doContinuum:
#            if (thisIon.Z - thisIon.Ion) in [0, 1]:
#                thisIon.twoPhoton(wavelength)
#                outList.append(thisIon.TwoPhoton)
#        outQueue.put(outList)
#    return
