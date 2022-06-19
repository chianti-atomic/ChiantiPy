import copy
from datetime import datetime

import numpy as np

import ChiantiPy
import ChiantiPy.tools.data as chdata
import ChiantiPy.tools.util as util
import ChiantiPy.tools.io as chio
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
    def __init__(self, temperature, eDensity, wvlRange, elementList=None, ionList=None, minAbund=None,
        keepIons=0, em=None, abundance=None, verbose=0, allLines=True):
        """
        Calculate the emission line spectrum as a function of temperature and density.

        """
        #
        t1 = datetime.now()
        # creates Intensity dict from first ion calculated
        setupIntensity = 0
        #
        self.argCheck(temperature=temperature, eDensity=eDensity, pDensity=None, em=em, verbose=verbose)
        self.Defaults=chdata.Defaults
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
        #
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
                        self.Intensity[akey] = np.hstack((copy.copy(self.Intensity[akey]), thisIon.Intensity[akey]))
                else:
                    setupIntensity = 1
#                                print(' creating Intensity dict from ion %s'%(ionS))
                    self.Intensity  = thisIon.Intensity
            else:
                if verbose:
                    print(thisIon.Intensity['errorMessage'])
        #
        #
        t2 = datetime.now()
        dt=t2-t1
        print(' elapsed seconds = %12.3f'%(dt.seconds))


