import copy
from datetime import datetime
import numpy as np
import ChiantiPy
import ChiantiPy.tools.data as chdata
import ChiantiPy.tools.constants as const
import ChiantiPy.tools.filters as chfilters
import ChiantiPy.tools.util as util
import ChiantiPy.tools.io as chio

from ChiantiPy.base import ionTrails
from ChiantiPy.base import specTrails


class spectrum(ionTrails, specTrails):
    '''
    Calculate the emission spectrum as a function of temperature and density.

    one of the convenient things is that all of the instantiated ion classes, determined
    through such keywords as 'elementList', 'ionList', and 'minAbund' are kept in a
    dictionary self.IonInstances where self.IonInstances['mg_7'] is the class instance of
    ChiantiPy.core.ion for 'mg_7'.  All its methods and attributes are available.

    includes elemental abundances and ionization equilibria

    the set of abundances, a file in $XUVTOP/abundance, can be set with the keyword argument 'abundanceName'

    temperature and density can be arrays but, unless the size of either is unity (1),
    the two must have the same size

    the returned spectrum will be convolved with a filter of the specified width on the
    specified wavelength array

    the default filter is gaussianR with a resolving power of 1000.  Other filters,
    such as gaussian, box and lorentz, are available in ChiantiPy.tools.filters.  When
    using the box filter, the width should equal the wavelength interval to keep the units
    of the continuum and line spectrum the same.

    Inherited methods include 'intensityList', 'intensityRatio' (between lines of different ions),
    'intensityRatioSave' and 'convolve'

    A selection of elements can be make with elementList a list containing the names of elements
    that are desired to be included, e.g., ['fe','ni']

    A selection of ions can be make with ionList containing the names of
    the desired lines in CHIANTI notation, i.e. C VI = c_6

    Both elementList and ionList can not be specified at the same time

    a minimum abundance can be specified so that the calculation can be speeded up
    by excluding elements with a low abundance. The default of minAbund is 1.e-6

    It is necessary to specify at least an elementList, an ionList, or a minAbund to select any ions
    for a spectrum calculation

    With solar photospheric abundances

    setting minAbund = 1.e-4 will include H, He, C, O, Ne

    setting minAbund = 2.e-5 adds  N, Mg, Si, S, Fe

    setting minAbund = 1.e-6 adds  Na, Al, Ar, Ca, Ni

    Setting doLines = 0 will skip the calculation of spectral lines.
    Setting doContinuum =0 will skip the continuum calculation.

    Setting em [for emission measure] will multiply the spectrum at each temperature
    by the value of em.

    em [for emission measure] can be a float or an array of the same length as the
    temperature/density

    keepIons  set this to keep the ion instances that have been calculated in a dictionary
    self.IonInstances with the keywords being the CHIANTI-style ion names

    abundance -  to select a particular set of abundances, set abundance to the name of a
    CHIANTI abundance file, without the '.abund' suffix, e.g. 'sun_photospheric_1998_grevesse'

    If set to a blank (''), a gui selection menu will popup and allow the selection of an
    set of abundances

    Parameters
    --------------

    temperature: `float`, `list`, `ndarray`
        the temperature(s) in K

    eDensity: float, ndarray
        eDensity: electron density in :math:`\mathrm{cm^{-3}}`

    wavelength:  `list` or `ndarray`
        wavelength:  array of  wavelengths, generally in Angstroms

    elementList:  `list`
        elementList:  list of elements to include, such as 'fe', 'ne', 's'

    ionList:  `list`
        ionList:  list of ions to include, such as 'fe_16', 'ne_10'

    minAbund:  `float`
        minAbund:  minimum abundance (relative to H) to include

    doLines:  `bool1
        doLines: if true, line intensities are calculated

    doContinuum:  `bool`
        doContinuum:  if true, continuum intensities are calculated only if wavelengths are in angstroms

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
    def __init__(self, temperature, eDensity, wavelength, filter=(chfilters.gaussianR, 1000.), label=None,
        elementList = None, ionList = None, minAbund=None, doLines=1, doContinuum=1, em=None, keepIons=0,
        abundance=None, verbose=0, allLines=1):
        #
        self.Defaults=chdata.Defaults

        if doContinuum and self.Defaults['wavelength'] != 'angstrom':
            print(' the continuum can only be calculated for wavelengths in angstroms')
            print(' set doContuum = False to continue')
            return

        wavelength = np.atleast_1d(wavelength)

        if wavelength.size < 2:
            print(' wavelength must have at least two values, current length %3i'%(wavelength.size))
            return

        t1 = datetime.now()
        # creates Intensity dict from first ion calculated
        setupIntensity = False
        #
        self.Wavelength = np.asarray(wavelength,  np.float64)
        self.WvlRange = np.asarray([self.Wavelength.min(),  self.Wavelength.max()],  np.float64)
        #
        self.argCheck(temperature=temperature, eDensity=eDensity, pDensity=None,  em=em,  verbose=verbose)

        nTempDens = self.NTempDens

        self.Labels = util.units(chdata.Defaults)

        xlabel = self.Labels['xlabel']
        ylabel = self.Labels['spectrumYlabel']

        if np.array_equal(self.Em, np.ones_like(self.Em)):
            ylabel += '($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'
        #
        if abundance is not None:
            ab = chio.abundanceRead(abundance)
            abundAll = ab['abundance']
            self.AbundanceName = abundance
        else:
            self.AbundanceName = self.Defaults['abundfile']
            abundAll = chdata.Abundance[self.AbundanceName]['abundance']

        # needed by ionGate
        self.AbundAll = abundAll
        self.Abundance = abundAll
        #
        self.MinAbund = minAbund
        wavelength = np.asarray(wavelength)
        nWvl = wavelength.size
        self.Wavelength = wavelength
        #
        freeFree = np.zeros((nTempDens, nWvl), np.float64).squeeze()
        freeBound = np.zeros((nTempDens, nWvl), np.float64).squeeze()
        twoPhoton = np.zeros((nTempDens, nWvl), np.float64).squeeze()
        lineSpectrum = np.zeros((nTempDens, nWvl), np.float64).squeeze()
        #
        self.IonsCalculated = []
        if keepIons:
            self.IonInstances = {}
            self.FfInstances = {}
            self.FbInstances = {}
        self.Finished = []
        #
        self.ionGate(elementList = elementList, ionList = ionList, minAbund=minAbund, doLines=doLines,
            doContinuum=doContinuum, verbose = verbose)
        #
        for akey in sorted(self.Todo.keys()):
            zStuff = util.convertName(akey)
            Z = zStuff['Z']
            ionstage = zStuff['Ion']
            dielectronic = zStuff['Dielectronic']
            abundance = self.Abundance[Z - 1]
            if verbose:
                print(' %5i %5s abundance = %10.2e '%(Z, const.El[Z-1],  abundance))
            if verbose:
                print(' doing ion %s for the following processes %s'%(akey, self.Todo[akey]))
            if 'ff' in self.Todo[akey]:
                if verbose:
                    print(' calculating ff continuum for :  %s'%(akey))
                FF = ChiantiPy.core.continuum(akey, temperature, abundance=abundance, em=em,  verbose=verbose)
                FF.freeFree(wavelength)
                freeFree += FF.FreeFree['intensity'].squeeze()
                if keepIons:
                    self.FfInstances[akey] = copy.deepcopy(FF)

            if 'fb' in self.Todo[akey]:
                if verbose:
                    print(' calculating fb continuum for :  %s'%(akey))
                FB = ChiantiPy.core.continuum(akey, temperature, abundance=abundance, em=em, verbose=verbose)
                FB.freeBound(wavelength)
                if 'errorMessage' not in FB.FreeBound.keys():
                    freeBound += FB.FreeBound['intensity'].squeeze()
                    if keepIons:
                        self.FbInstances[akey] = copy.deepcopy(FB)
                else:
                    if verbose:
                        print(FB.FreeBound['errorMessage'])
            if 'line' in self.Todo[akey]:
                if verbose:
                    print(' calculating spectrum for  :  %s'%(akey))
                thisIon = ChiantiPy.core.ion(akey, temperature, eDensity, pDensity='default', abundance=abundance, em=em, verbose=verbose)
                thisIon.intensity(allLines=allLines)
                self.IonsCalculated.append(akey)
                if 'errorMessage' not in  list(thisIon.Intensity.keys()):
                    self.Finished.append(akey)
                    thisIon.spectrum(wavelength, filter=filter, allLines=allLines)
                    if keepIons:
                        self.IonInstances[akey] = copy.deepcopy(thisIon)
                    if setupIntensity:
                        for bkey in self.Intensity:
                            self.Intensity[bkey] = np.hstack((copy.copy(self.Intensity[bkey]),
                                thisIon.Intensity[bkey]))
                    else:
                        setupIntensity = True
                        self.Intensity  = thisIon.Intensity
                    lineSpectrum += thisIon.Spectrum['intensity'].squeeze()
                else:
                    if verbose:
                        print(thisIon.Intensity['errorMessage'])
                # get 2 photon emission for H and He sequences
                if (Z - ionstage) in [0, 1] and not dielectronic:
                    thisIon.twoPhoton(wavelength)
                    twoPhoton += thisIon.TwoPhoton['intensity'].squeeze()

        self.FreeFree = {'wavelength':wavelength, 'intensity':freeFree.squeeze()}
        self.FreeBound = {'wavelength':wavelength, 'intensity':freeBound.squeeze()}
        self.LineSpectrum = {'wavelength':wavelength, 'intensity':lineSpectrum.squeeze()}
        self.TwoPhoton = {'wavelength':wavelength, 'intensity':twoPhoton.squeeze()}
        cont = freeFree.squeeze() + freeBound.squeeze() + twoPhoton.squeeze()
        self.Continuum = {'wavelength':wavelength, 'intensity':cont}
        #
        #
        total = freeFree + freeBound + lineSpectrum + twoPhoton
        self.Total = total

        t2 = datetime.now()
        dt=t2-t1
        print(' elapsed seconds = %12.3f'%(dt.seconds))
        if nTempDens == 1:
            integrated = total
        else:
            integrated = total.sum(axis=0)
        #
        if type(label) == type(''):
            if hasattr(self, 'Spectrum'):
                self.Spectrum[label] = {'wavelength':wavelength, 'intensity':total.squeeze(),
                    'filter':filter[0],   'filterWidth':filter[1], 'integrated':integrated, 'em':em,
                    'ions':self.IonsCalculated, 'Abundance':self.AbundanceName, 'xlabel':xlabel,
                    'ylabel':ylabel, 'minAbund':minAbund}

            else:
                self.Spectrum = {label:{'wavelength':wavelength, 'intensity':total.squeeze(),
                    'filter':filter[0],   'filterWidth':filter[1], 'integrated':integrated, 'em':em,
                    'ions':self.IonsCalculated, 'Abundance':self.AbundanceName, 'xlabel':xlabel,
                    'ylabel':ylabel, 'minAbund':minAbund}}

        else:
            self.Spectrum = {'wavelength':wavelength, 'intensity':total.squeeze(),
                'filter':filter[0],   'filterWidth':filter[1], 'integrated':integrated,
                'ions':self.IonsCalculated, 'Abundance':self.AbundanceName, 'xlabel':xlabel,
                'ylabel':ylabel, 'minAbund':minAbund}




