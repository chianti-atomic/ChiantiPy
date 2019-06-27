from datetime import datetime
import copy
import warnings

import numpy as np
try:
    from ipyparallel import Client
except ImportError:
    warnings.warn("ipyparallel not found. You won't be able to use the ipymspectrum module")

import ChiantiPy
import ChiantiPy.tools.data as chdata
import ChiantiPy.tools.constants as const
import ChiantiPy.tools.filters as chfilters
import ChiantiPy.tools.util as util
import ChiantiPy.Gui as chGui
from ChiantiPy.base import ionTrails
from ChiantiPy.base import specTrails


class ipymspectrum(ionTrails, specTrails):
    '''
    this is the multiprocessing version of spectrum for using inside an IPython Qtconsole or notebook.

    be for creating an instance, it is necessary to type something like the following into a console

    > ipcluster start   --n=3

    this is the way to invoke things under the IPython 6 notation

    Calculate the emission spectrum as a function of temperature and density.

    temperature and density can be arrays but, unless the size of either is one (1),
    the two must have the same size

    the returned spectrum will be convolved with a filter of the specified width on the
    specified wavelength array

    the default filter is gaussianR with a resolving power of 100.  Other filters,
    such as gaussian, box and lorentz, are available in ChiantiPy.filters.  When using the box filter,
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
    def __init__(self, temperature, eDensity, wavelength, filter=(chfilters.gaussianR, 1000.), label=None, elementList = None, ionList = None, minAbund=None, keepIons=0, doLines=1, doContinuum=1, allLines = 1, em=None, abundance=None, verbose=0,  timeout=0.1):
        #
        wavelength = np.atleast_1d(wavelength)
        if wavelength.size < 2:
            print(' wavelength must have at least two values, current length %3i'%(wavelength.size))
            return

        t1 = datetime.now()
        #
        rcAll = Client()
#        all_engines = rcAll[:]
        lbvAll = rcAll.load_balanced_view()
        #
        #
        # creates Intensity dict from first ion calculated
        #
        setupIntensity = 0
        #
        self.Defaults = chdata.Defaults
        #
        self.argCheck(temperature=temperature, eDensity=eDensity, pDensity = None,  em=em)
        if verbose:
            print('NTempDens:  %5i'%(self.NTempDens))
            #
        #
        if self.Em.max() == 1.:
            ylabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1} \AA^{-1}$ ($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'
        else:
            ylabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1} \AA^{-1}$ $'
        if verbose:
            print('len of self.Em %5i'%(len(self.Em)))
        #
        #
        #
        xlabel = 'Wavelength ('+self.Defaults['wavelength'] +')'
        #
        self.AllLines = allLines
        #
        #
        if abundance != None:
            if abundance in list(chdata.Abundance.keys()):
                self.AbundanceName = abundance
            else:
                abundChoices = list(chdata.Abundance.keys())
                abundChoice = chGui.gui.selectorDialog(abundChoices,label='Select Abundance name', multiChoice=False)
                abundChoice_idx = abundChoice.selectedIndex
                self.AbundanceName = abundChoices[abundChoice_idx[0]]
                print((' Abundance chosen:  %s '%(self.AbundanceName)))
        else:
            self.AbundanceName = self.Defaults['abundfile']
        #
        #
        abundAll = chdata.Abundance[self.AbundanceName]['abundance']
        self.AbundAll = abundAll
        self.MinAbund = minAbund
        #
        #ionInfo = chio.masterListInfo()
        wavelength = np.asarray(wavelength)
        nWvl = wavelength.size
        self.Wavelength = wavelength
        ntemp = self.Ntemp
        #
        #
        freeFree = np.zeros((ntemp, nWvl), np.float64).squeeze()
        freeBound = np.zeros((ntemp, nWvl), np.float64).squeeze()
        twoPhoton = np.zeros((self.NTempDens, nWvl), np.float64).squeeze()
        lineSpectrum = np.zeros((self.NTempDens, nWvl), np.float64).squeeze()
        #
         #
        allInpt = []
        #
        if keepIons:
            self.IonInstances = {}
            self.FbInstances = {}
            self.FfInstances = {}
        #
        # ionGate creates the self.Todo list
        #
        self.ionGate(elementList = elementList, ionList = ionList, minAbund=minAbund, doLines=doLines, doContinuum=doContinuum, verbose = verbose)
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
                allInpt.append([akey, 'ff', temperature, wavelength, abundance, em])
            if 'fb' in self.Todo[akey]:
                allInpt.append([akey, 'fb', temperature, wavelength, abundance, em])
            if 'line' in self.Todo[akey]:
                allInpt.append([akey, 'line', temperature, eDensity, wavelength, filter, allLines, abundance, em, doContinuum])
        #
        result = lbvAll.map_sync(doAll, allInpt)
        if verbose:
            print(' got all ff, fb, line results')
        ionsCalculated = []
        #
        for ijk in range(len(result)):
            out = result[ijk]
            if type(out) != list:
                print(' a problem has occured - this can be caused by')
                print('running Python3 and not using ipcluster3')
                return
            ionS = out[0]
            if verbose:
                print(' collecting calculation for %s'%(ionS))
            ionsCalculated.append(ionS)
            calcType = out[1]
            if verbose:
                print(' processing %s results'%(calcType))
            #
            if calcType == 'ff':
                thisFf = out[2]
                if keepIons:
                    self.FfInstances[ionS] = thisFf
                freeFree += thisFf['intensity']
            elif calcType == 'fb':
                thisFb = out[2]
                if verbose:
                    print(' fb ion = %s'%(ionS))
                if 'intensity' in thisFb.keys():
                    if 'errorMessage' not in sorted(thisFb.keys()):
                        if keepIons:
                            self.FbInstances[ionS] = thisFb
                        freeBound += thisFb['intensity']
                    else:
                        print(thisFb['errorMessage'])
            elif calcType == 'line':
                thisIon = out[2]
                if not 'errorMessage' in sorted(thisIon.Intensity.keys()):
                    if keepIons:
                        self.IonInstances[ionS] = thisIon
                    thisIntensity = thisIon.Intensity
    ##                self.IonInstances.append(copy.deepcopy(thisIon))
                    if setupIntensity:
                        for akey in sorted(self.Intensity.keys()):
                            self.Intensity[akey] = np.hstack((self.Intensity[akey], thisIntensity[akey]))
                    else:
                        setupIntensity = 1
                        self.Intensity  = thisIntensity
                    #
                    lineSpectrum += thisIon.Spectrum['intensity']
                   # check for two-photon emission
                    if len(out) == 4:
                        tp = out[3]
                        if verbose:
                            for akey in tp:
                                print(' tp keys %s'%(akey))
                        if self.NTempDens == 1:
                            twoPhoton += tp['intensity'].squeeze()
                        else:
                            for iTempDen in range(self.NTempDens):
                                twoPhoton[iTempDen] += tp['intensity'][iTempDen]
                else:
                    if 'errorMessage' in sorted(thisIon.Intensity.keys()):
                        print(thisIon.Intensity['errorMessage'])
        #
        #
        self.IonsCalculated = ionsCalculated
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
        print(' elapsed seconds = %12.3e'%(dt.seconds))
        rcAll.purge_results('all')
        #
        if self.NTempDens == 1:
            integrated = total
        else:
            integrated = total.sum(axis=0)
        #
        if type(label) == type(''):
            if hasattr(self, 'Spectrum'):
                print(' hasattr = true')
                self.Spectrum[label] = {'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':em,  'Abundance':self.AbundanceName,
                                            'xlabel':xlabel, 'ylabel':ylabel}
            else:
                self.Spectrum = {label:{'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':em, 'Abundance':self.AbundanceName,
                                'xlabel':xlabel, 'ylabel':ylabel}}
        else:
            self.Spectrum = {'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':em, 'Abundance':self.AbundanceName,
                                'xlabel':xlabel, 'ylabel':ylabel}
    #
    # -------------------------------------------------------------------------
    #
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
        em = inpt[5]
        FF = ChiantiPy.core.continuum(ionS, temperature, abundance=abund, em=em)
        FF.freeFree(wavelength)
        # can not do a deep copy of
#        return [ionS, calcType, copy.deepcopy(cont)]
        return [ionS, calcType, copy.copy(FF.FreeFree)]
    elif calcType == 'fb':
        temperature = inpt[2]
        wavelength = inpt[3]
        abund = inpt[4]
        em = inpt[5]
        cont = ChiantiPy.core.continuum(ionS, temperature, abundance=abund, em=em)
        #try:
        cont.freeBound(wavelength)
            #return [ionS, calcType, {'rate': cont.FreeBound}]
        return [ionS, calcType, copy.copy(cont.FreeBound)]
        #except ValueError:
            #return [ionS, calcType, {'errorMessage': 'No free-bound information available.'}]
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
        thisIon = ChiantiPy.core.ion(ionS, temperature, density, abundance=abund, em=em)
        thisIon.intensity(allLines = allLines)
        if 'errorMessage' not in thisIon.Intensity.keys():
            thisIon.spectrum(wavelength,  filter=filter, allLines=allLines)
        outList = [ionS, calcType, copy.deepcopy(thisIon)]
        if not thisIon.Dielectronic and doContinuum:
            if (thisIon.Z - thisIon.Ion) in [0, 1]:
                thisIon.twoPhoton(wavelength)
                outList.append(thisIon.TwoPhoton)
        return outList
