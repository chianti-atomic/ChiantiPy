from datetime import datetime
import copy
import multiprocessing as mp

import numpy as np

import ChiantiPy.tools.data as chdata
import ChiantiPy.tools.constants as const
import ChiantiPy.tools.filters as chfilters
import ChiantiPy.tools.util as util
import ChiantiPy.Gui as chGui
from ChiantiPy.base import ionTrails
from ChiantiPy.base import specTrails
import ChiantiPy.tools.mputil as mputil


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
    def __init__(self, temperature, eDensity, wavelength, filter=(chfilters.gaussianR, 1000.), label=0, elementList = None, ionList = None, minAbund=None, keepIons=0, abundance=None,  doLines=1, doContinuum=1, allLines = 1, em=None,  proc=3, verbose = 0,  timeout=0.1):
        #
        wavelength = np.atleast_1d(wavelength)
        if wavelength.size < 2:
            print(' wavelength must have at least two values, current length %3i'%(wavelength.size))
            return
        t1 = datetime.now()
        # creates Intensity dict from first ion calculated
        setupIntensity = 0
        #
        self.Defaults = chdata.Defaults
        #

        self.argCheck(temperature=temperature, eDensity=eDensity, pDensity=None, em=em)

        nTempDens = self.NTempDens

        if self.Em.max() == 1.:
            ylabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1} \AA^{-1}$ ($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'
        else:
            ylabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1} \AA^{-1}$ $'
        #
        #
        if self.Defaults['wavelength'] == 'angstrom':
            xlabel = 'Wavelength ('+self.Defaults['wavelength'].capitalize() +')'
        else:
            xlabel = 'Wavelength ('+self.Defaults['wavelength'] +')'
        #
        #
        self.AllLines = allLines
        self.MinAbund = minAbund
        #
        if abundance is not None:
            if type(abundance) == str:
                if abundance in chdata.AbundanceList:
                    self.AbundanceName = abundance
                else:
                    abundChoices = chdata.AbundanceList
                    abundChoice = chGui.gui.selectorDialog(abundChoices,label='Select Abundance name', multiChoice=False)
                    abundChoice_idx = abundChoice.selectedIndex
                    self.AbundanceName = abundChoices[abundChoice_idx[0]]
                    if verbose:
                        print('Abundance file chosen: %s'%(self.AbundanceName))
            else:
                print(' keyword abundance must be a string, either a blank (\'\') or the name of an abundance file')
                return
        else:
            self.AbundanceName = self.Defaults['abundfile']
        if hasattr(self,'AbundanceName'):
            self.Abundance = chdata.Abundance[self.AbundanceName]['abundance']
        #
        abundAll = chdata.Abundance[self.AbundanceName]['abundance']
        self.AbundAll = abundAll
        #
        wavelength = np.asarray(wavelength)
        nWvl = wavelength.size
        self.Wavelength = wavelength
        #
        proc = min([proc, mp.cpu_count()])
        #
        freeFree = np.zeros((nTempDens, nWvl), np.float64).squeeze()
        freeBound = np.zeros((nTempDens, nWvl), np.float64).squeeze()
        twoPhoton = np.zeros((nTempDens, nWvl), np.float64).squeeze()
        lineSpectrum = np.zeros((nTempDens, nWvl), np.float64).squeeze()
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
            if doLines:
                self.IonInstances = {}
        self.Finished = []
        #

        self.ionGate(elementList = elementList, ionList = ionList, minAbund=minAbund, doLines=doLines, doContinuum=doContinuum, verbose = verbose)
        for one in self.Todo.keys():
            print(' %s  %s'%(one, self.Todo[one]))
        #
        for akey in sorted(self.Todo.keys()):
#            zStuff = util.convertName(akey)
#            Z = zStuff['Z']
#            abundance = self.Abundance[Z - 1]
#            if verbose:
#                print(' %5i %5s abundance = %10.2e '%(Z, const.El[Z-1],  abundance))
            if verbose:
                print(' doing ion %s for the following processes %s'%(akey, self.Todo[akey]))
            if 'ff' in self.Todo[akey]:
                ffWorkerQ.put((akey, temperature, wavelength, abundance, em))
            if 'fb' in self.Todo[akey]:
                fbWorkerQ.put((akey, temperature, wavelength, abundance, em))
            if 'line' in self.Todo[akey]:
                ionWorkerQ.put((akey, temperature, eDensity, wavelength, filter, allLines, abundance, em, doContinuum))
        #
        ffWorkerQSize = ffWorkerQ.qsize()
        fbWorkerQSize = fbWorkerQ.qsize()
        ionWorkerQSize = ionWorkerQ.qsize()
        if doContinuum:
            ffProcesses = []
            for i in range(proc):
                p = mp.Process(target=mputil.doFfQ, args=(ffWorkerQ, ffDoneQ))
                p.start()
                ffProcesses.append(p)
    #       timeout is not necessary
            for p in ffProcesses:
                if p.is_alive():
                    p.join(timeout=timeout)
            #
            for iff in range(ffWorkerQSize):
                thisFreeFree = ffDoneQ.get()
                freeFree += thisFreeFree['intensity']
            for p in ffProcesses:
                if not isinstance(p, str):
                    p.terminate()
        #
            fbProcesses = []
            for i in range(proc):
                p = mp.Process(target=mputil.doFbQ, args=(fbWorkerQ, fbDoneQ))
                p.start()
                fbProcesses.append(p)
    #       timeout is not necessary
            for p in fbProcesses:
                if p.is_alive():
                    p.join(timeout=timeout)
            #
            for ifb in range(fbWorkerQSize):
                thisFreeBound = fbDoneQ.get()
                freeBound += thisFreeBound['intensity'].squeeze()

            for p in fbProcesses:
                if not isinstance(p, str):
                    p.terminate()
        #
        if doLines:
            ionProcesses = []
            if ionWorkerQSize < proc:
                proc = ionWorkerQSize
            for i in range(proc):
                p = mp.Process(target=mputil.doIonQ, args=(ionWorkerQ, ionDoneQ))
                p.start()
                ionProcesses.append(p)
        #       timeout is not necessary
            for p in ionProcesses:
        #            if p.is_alive():
        #                p.join()
                    p.join(timeout=timeout)
            #
            for ijk in range(ionWorkerQSize):
                out = ionDoneQ.get()
                ionS = out[0]
                if verbose:
                    print(' collecting ion calculation for %s'%(ionS))
                thisIon = out[1]
                thisIntensity = thisIon.Intensity
                if not 'errorMessage' in sorted(thisIntensity.keys()):
                    self.Finished.append(ionS)
                    if keepIons:
                        self.IonInstances[ionS] = copy.deepcopy(thisIon)
                    if setupIntensity:
                        for akey in sorted(self.Intensity.keys()):
                            self.Intensity[akey] = np.hstack((copy.copy(self.Intensity[akey]), thisIntensity[akey]))
                    else:
                        setupIntensity = 1
                        self.Intensity  = thisIntensity
                    #
                    if not 'errorMessage' in sorted(thisIon.Spectrum.keys()):
                        lineSpectrum += thisIon.Spectrum['intensity']
                   # check for two-photon emission
                    if len(out) == 3:
                        tp = out[2]
                        twoPhoton += tp['intensity'].squeeze()
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
        if nTempDens == 1:
            integrated = total
        else:
            integrated = total.sum(axis=0)
        #
        if type(label) == type(''):
            if hasattr(self, 'Spectrum'):
                self.Spectrum[label] = {'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'ions':self.IonsCalculated, 'em':em,
                'Abundance':self.AbundanceName, 'xlabel':xlabel, 'ylabel':ylabel, 'minAbund':minAbund}
            else:
                self.Spectrum = {label:{'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'ions':self.IonsCalculated, 'em':em,
                'Abundance':self.AbundanceName, 'xlabel':xlabel, 'ylabel':ylabel}, 'minAbund':minAbund}
        else:
            self.Spectrum ={'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'ions':self.IonsCalculated,
            'em':em, 'Abundance':self.AbundanceName, 'xlabel':xlabel, 'ylabel':ylabel, 'minAbund':minAbund}
