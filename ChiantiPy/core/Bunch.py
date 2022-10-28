import os.path
import copy
from datetime import datetime

import numpy as np

import ChiantiPy
import ChiantiPy.tools.data as chdata
import ChiantiPy.tools.util as util
import ChiantiPy.tools.io as chio
import ChiantiPy.tools.filters as chfilters
import ChiantiPy.Gui as chGui
from ChiantiPy.base import ionTrails
from ChiantiPy.base import specTrails


class bunch(ionTrails, specTrails):
    '''
    Calculate the emission line spectrum as a function of temperature and density.

    'bunch' is very similar to 'spectrum' except that continuum is not calculated and
    the spectrum is not convolved over a filter.  However, this can be done with the
    inherited convolve method

    one of the convenient things is that all of the instantiated ion classes,
    determined through such keywords as 'elementList', 'ionList', and 'minAbund' are
    kept in a dictionary self.IonInstances where self.IonInstances['mg_7'] is the
    class instance of ChiantiPy.core.ion for 'mg_7'.  All its methods and attributes
    are available.

    includes elemental abundances and ionization equilibria

    the set of abundances, a file in $XUVTOP/abundance, can be set with the keyword
    argument 'abundanceName'

    temperature and density can be arrays but, unless the size of either is one (1),
    the two must have the same size

    Inherited methods include 'intensityList', 'intensityRatio' (between lines of different
    ions), and 'intensityRatioSave' and 'convolve'.

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

    At least one of elementList, ionList, or minAbund must be set in order for 'bunch' to include
    any ions.

    Setting em will multiply the spectrum at each temperature by the value of em.

    em [for emission measure], can be a float or an array of the same length as the
    temperature/density

    Parameters
    ----------

    temperature: `float`, `list`, `ndarray`
        the temperature(s) in K

    eDensity: float, ndarray
        eDensity: electron density in :math:`\mathrm{cm^{-3}}`

    wvlRange:  2 element `list` or `ndarray`
        wvlRange:  range of wavelengths to consider, generally in angstroms

    elementList:  `list`
        elementList:  list of elements to include, such as 'fe', 'ne', 's'

    ionList:  `list`
        ionList:  list of ions to include, such as 'fe_16', 'ne_10'

    minAbund:  `float`
        minAbund:  minimum abundance (relative to H) to include

    keepIons:  `bool`
        keepIons:  keep the ion instances used in the calculation
            should be used with caution otherwise the bunch instance
            can become quite large

    em:  `float`, `list`, `ndarray`
        em:  the emission measure

    abundance:  `str`
        abuncance:  the file name of the abuncance set to be used
            must be one in the $XUVTOP/abund directory

    allLInes:  `bool`
        allLines:  whether or not to include unobserved lines

    verbose:  `bool`
        verbose:  whether to allow certain print statements

    '''
    #
    # ------------------------------------------------------------------------------------
    #
    def __init__(self, temperature, eDensity, wvlRange, elementList=None, ionList=None,
        minAbund=None, keepIons=False, em=None, abundance=None, verbose=False, allLines=True):
        """
        Calculate the emission line spectrum as a function of temperature and density.

        """

        #
        t1 = datetime.now()
        # creates Intensity dict from first ion calculated
        setupIntensity = False
        #
        self.argCheck(temperature=temperature, eDensity=eDensity, pDensity=None, em=em, verbose=verbose)
        self.Defaults=chdata.Defaults
        self.Labels = util.units(chdata.Defaults)
        #
        #
        if abundance is not None:
            ab = chio.abundanceRead(abundance)
            abundAll = ab['abundance']
            self.AbundanceName = ab['abundancename']
        else:
            self.AbundanceName = self.Defaults['abundfile']
        #
            abundAll = chdata.Abundance[self.AbundanceName]['abundance']
        # needed by ionGate
        self.AbundAll = abundAll
        self.Abundance = abundAll
        #
#        # this is usually done in argCheck
#        if em is not None:
#            em = np.atleast_1d(em)
#            self.Em = em
#            if em.size == 1:
#                self.Em = np.tile(em,self.NTempDens)
#
#            elif em.size != self.NTempDens:
#                raise ValueError('the size of em must be either 1 or the size of the larger of temperature or density %5i'%(self.NTempDens))
#        else:
#            self.Em = np.ones_like(self.Temperature, np.float64)

#        if self.Defaults['wavelength'] == 'angstrom':
#            xlabel = 'Wavelength \u212B'
#        else:
#            xlabel = 'Wavelength ('+self.Defaults['wavelength'] +')'
#
#        # unicode character for angstrom is \u212B
#        if self.Em.max() == 1.:
#            ylabel = 'erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$ ($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'
#        else:
#            ylabel = 'erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$'

        xlabel = self.Labels['xlabel']
        ylabel = self.Labels['intensityYlabel']

        if np.array_equal(self.Em, np.ones_like(self.Em)):
            ylabel += '($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'

#        nonzed = abundAll > 0.
#        minAbundAll = abundAll[nonzed].min()
#        # if minAbund is even set
#        if minAbund:
#            if minAbund < minAbundAll:
#                minAbund = minAbundAll
#        self.minAbund = minAbund
#        ionInfo = chio.masterListInfo()
        #        #
        self.IonsCalculated = []
        if keepIons:
            self.IonInstances = {}
        self.Finished = []
        #
        # also needed by ionGate
        self.WvlRange = np.asarray(wvlRange, 'float64')
        #
        self.ionGate(elementList = elementList, ionList = ionList, minAbund=minAbund, doLines=1,
            doContinuum=0, verbose = False)
        #
        for ionS in sorted(self.Todo.keys()):
            nameStuff = util.convertName(ionS)
            Z = nameStuff['Z']

            if verbose:
                print(' calculating %s'%(ionS))
            thisIon = ChiantiPy.core.ion(ionS, temperature, eDensity, abundance=abundAll[Z-1],  em=em)
            thisIon.intensity(allLines = allLines)
            self.IonsCalculated.append(ionS)
            #
            if 'errorMessage' not in  list(thisIon.Intensity.keys()):
                self.Finished.append(ionS)
#                            thisIon.spectrum(wavelength, filter=filter)
                if keepIons:
                    self.IonInstances[ionS] = copy.deepcopy(thisIon)
                if setupIntensity:
                    for akey in self.Intensity:
                        self.Intensity[akey] = np.hstack((copy.copy(self.Intensity[akey]),
                            thisIon.Intensity[akey]))
                else:
                    setupIntensity = True
#                                print(' creating Intensity dict from ion %s'%(ionS))
                    self.Intensity  = thisIon.Intensity
            else:
                if verbose:
                    print(thisIon.Intensity['errorMessage'])
        #
        self.Intensity['xlabel'] = xlabel
        self.Intensity['ylabel'] = ylabel
        #
        t2 = datetime.now()
        dt=t2-t1
        print(' elapsed seconds = %12.3f'%(dt.seconds))


    def convolve(self, wavelength=None, filter=(chfilters.gaussianR, 1000.), label=None, verbose=False):
        '''
        the first application of spectrum calculates the line intensities within the specified wavelength
        range and for set of ions specified

        wavelength will not be used if applied to 'spectrum' objects

        wavelength IS need for 'bunch' objects - in this case, the wavelength should not extend beyond
        the limits of the wvlRange used for the 'bunch' calculation

        Keyword Arguments
        -----------------

        wavelength:  'int', `list`
            if an `int`, the attribute 'Wavelength' is looked for
            otherwise, wavelength is used

        filter: `tuple`
            first elements if one of the ChiantiPy.tools.filters object
            second element is the width appropriate to the filter

        label:  `str`
            if set, creates a Spectrum[label] attribute

        verbose: `bool`
            if True, prints info to the terminal

        '''
        useFilter = filter[0]
        useFactor = filter[1]

        #
        if label is not None:
            if not isinstance(label, str):
                print(' label must either be None or a string')
        #
        t1 = datetime.now()

        # unicode character for angstrom is \u212B
        if self.Em.max() == 1.:
            ylabel = self.Labels['spectrumYlabel'] + ' ($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'
        else:
            ylabel = self.Labels['spectrumYlabel']

        xlabel = self.Labels['xlabel']

        #:
        if hasattr(self, 'Wavelength'):
                wavelength = self.Wavelength
        elif wavelength is not None:
            self.Wavelength = wavelength
        else:
            print(' a wavelength array must be given')
            return

        if not hasattr(self, 'NTempDens'):
            self.NTempDens = max([self.Ntemp,  self.Ndens])

#        aspectrum = np.zeros((self.NTempDens, wavelength.size), np.float64)

        if hasattr(self, 'IonInstances'):
            for akey in sorted(self.IonInstances.keys()):
                if verbose:
                    print( ' trying ion = %s'%(akey))
                if not 'errorMessage' in sorted(self.IonInstances[akey].Intensity.keys()):
                    if verbose:
                        print(' doing convolve on ion %s '%(akey))
                    self.IonInstances[akey].spectrum(wavelength, filter)
                    if verbose:
                        if 'errorMessage' in sorted(self.IonInstances[akey].Spectrum.keys()):
                            print(self.IonInstances[akey].Spectrum['errorMessage'])
#                    else:
#                        aspectrum += self.IonInstances[akey].Spectrum['intensity']
    #                if self.NTempDens == 1:
    #                    lineSpectrum += thisIon.Spectrum['intensity']
    #                else:
    #                    for iTempDen in range(self.NTempDens):
    #                        lineSpectrum[iTempDen] += thisIon.Spectrum['intensity'][iTempDen]
                else:
                    if verbose:
                        if 'errorMessage' in sorted(self.IonInstances[akey].Intensity.keys()):
                            print(self.IonInstances[akey].Intensity['errorMessage'])

        bspectrum = np.zeros((self.NTempDens, wavelength.size), np.float64)
        if hasattr(self, 'Intensity'):
            for itemp in range(self.NTempDens):
                for iwvl, awvl in enumerate(self.Intensity['wvl']):
                    bspectrum[itemp] += useFilter(wavelength, awvl, factor=useFactor)*self.Intensity['intensity'][itemp, iwvl]

#        self.LineSpectrum = {'wavelength':wavelength, 'intensity':aspectrum.squeeze()}
        #
#        total = self.LineSpectrum['intensity']
        #
#        # the following is required in order to be applied to both a 'spectrum' and a 'bunch' object
#        #
#        if hasattr(self, 'FreeFree'):
#            total += self.FreeFree['intensity']
#        if hasattr(self, 'FreeBound'):
#            total += self.FreeBound['intensity']
#        if hasattr(self, 'TwoPhoton'):
#            total += self.TwoPhoton['intensity']
#        self.Total = total
        #
        #
        if self.NTempDens == 1:
            integrated = bspectrum
        else:
            integrated = bspectrum.sum(axis=0)
        #
        t2 = datetime.now()
        dt = t2 - t1
        print(' elapsed seconds = %12.3e'%(dt.seconds))
        #
        self.Spectrum = {'wavelength':wavelength, 'intensity':bspectrum.squeeze(),
            'integrated':integrated,'filter':filter[0], 'filterWidth':filter[1],
            'Abundance':self.AbundanceName, 'xlabel':xlabel, 'ylabel':ylabel}
#        if label is not None:
#        if hasattr(self, 'IonInstances'):
#            print(' has IonInstances')
#            if hasattr(self, 'Spectrum'):
#                self.Spectrum[label] = {'wavelength':wavelength, 'intensity':aspectrum.squeeze(),
#                    'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated,
#                    'em':self.Em,  'Abundance':self.AbundanceName, 'xlabel':xlabel, 'ylabel':ylabel}
#            else:
#                self.Spectrum = {label:{'wavelength':wavelength, 'intensity':aspectrum.squeeze(),
#                    'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated,
#                    'em':self.Em,  'Abundance':self.AbundanceName}, 'xlabel':xlabel, 'ylabel':ylabel}
#        else:


