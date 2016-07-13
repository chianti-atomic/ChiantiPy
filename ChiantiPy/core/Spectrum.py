import copy
from datetime import datetime
import numpy as np
import ChiantiPy
import ChiantiPy.tools.data as chdata
import ChiantiPy.tools.constants as const
import ChiantiPy.tools.filters as chfilters
import ChiantiPy.tools.util as util
import ChiantiPy.tools.io as chio
#
import ChiantiPy.Gui as chGui
#
from ._IonTrails import _ionTrails
from ._SpecTrails import _specTrails

defaults = chdata.Defaults

class spectrum(_ionTrails, _specTrails):
    '''
    Calculate the emission spectrum as a function of temperature and density.

    one of the convenient things is that all of the instantiated ion classes, determined through such keywords as 'elementList',
    'ionList', and 'minAbund' are kept in a dictionary self.IonInstances where self.IonInstances['mg_7'] is the class instance of
    chianti.core.ion for 'mg_7'.  All its methods and attributes are available.

    includes elemental abundances and ionization equilibria

    the set of abundances, a file in $XUVTOP/abundance, can be set with the keyword argument 'abundanceName'

    temperature and density can be arrays but, unless the size of either is one (1),
    the two must have the same size

    the returned spectrum will be convolved with a filter of the specified width on the
    specified wavelength array

    the default filter is gaussianR with a resolving power of 1000.  Other filters,
    such as gaussian, box and lorentz, are available in chianti.filters.  When using the box filter,
    the width should equal the wavelength interval to keep the units of the continuum and line
    spectrum the same.

    Inherited methods include 'intensityList', 'intensityRatio' (between lines of different ions), 'intensityRatioSave'
    and 'convolve'

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

    Setting doContinuum =0 will skip the continuum calculation.

    Setting em will multiply the spectrum at each temperature by the value of em.

    em [for emission measure], can be a float or an array of the same length as the
    temperature/density
    '''
    def __init__(self, temperature, eDensity, wavelength, filter=(chfilters.gaussianR, 1000.), label=0, elementList = 0, ionList = 0, minAbund=0, doContinuum=1, em=0, keepIons=0,  abundanceName=0, verbose=0, allLines=1):
        #
        t1 = datetime.now()
        # creates Intensity dict from first ion calculated
        setupIntensity = 0
        #
        masterlist = chdata.MasterList
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
        self.Defaults=defaults
        self.Temperature = np.asarray(temperature, 'float64')
        nTemp = self.Temperature.size
        self.EDensity = np.asarray(eDensity, 'float64')
        nDen = self.EDensity.size
        self.NTempDen = max([nTemp, nDen])
        nTempDen = self.NTempDen
        self.Wavelength = wavelength
        #
        if type(em) == int and em == 0:
            em = np.ones(self.NTempDen, 'float64')
        elif type(em) == float and em > 0.:
            em = np.ones(self.NTempDen, 'float64')*em
        elif type(em) == list or type(em) == tuple:
            em = np.asarray(em, 'float64')
        self.Em = em
        #
        if self.Em.any() > 0.:
            ylabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1} \AA^{-1}$ $'
        else:
            ylabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1} \AA^{-1}$ ($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'
        #
        xlabel = 'Wavelength ('+self.Defaults['wavelength'] +')'
        #
        #self.AbundanceName = defaults['abundfile']
        #self.AbundanceAll = chdata.AbundanceAll
        #
        if abundanceName:
            if abundanceName in list(chdata.Abundance.keys()):
                self.AbundanceName = abundanceName
            else:
                abundChoices = list(chdata.Abundance.keys())
#                for one in wvl[topLines]:
#                    wvlChoices.append('%12.3f'%(one))
                abundChoice = chGui.gui.selectorDialog(abundChoices,label='Select Abundance name')
                abundChoice_idx = abundChoice.selectedIndex
                self.AbundanceName = abundChoices[abundChoice_idx[0]]
                abundanceName = self.AbundanceName
                print((' Abundance chosen:  %s '%(self.AbundanceName)))
        else:
            self.AbundanceName = self.Defaults['abundfile']
        #
        abundAll = chdata.Abundance[self.AbundanceName]['abundance']
        # needed by ionGate
        self.AbundAll = abundAll
        #
        nonzed = abundAll > 0.
        minAbundAll = abundAll[nonzed].min()
        # if minAbund is even set
        if minAbund:
            if minAbund < minAbundAll:
                minAbund = minAbundAll
        else:
            minAbund = minAbundAll
        self.MinAbund = minAbund
#        ionInfo = chio.masterListInfo()
        wavelength = np.asarray(wavelength)
        nWvl = wavelength.size
        self.Wavelength = wavelength
        wvlRange = [wavelength.min(), wavelength.max()]
        #
        freeFree = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        freeBound = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        twoPhoton = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        lineSpectrum = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        #
        self.IonsCalculated = []
        if keepIons:
            self.IonInstances = {}
        self.Finished = []
        #
        self.ionGate(elementList = elementList, ionList = ionList, minAbund=minAbund, doContinuum=doContinuum, verbose = verbose)
        #
        for akey in sorted(self.Todo.keys()):
            zStuff = util.convertName(akey)
            Z = zStuff['Z']
            ionstage = zStuff['Ion']
            dielectronic = zStuff['Dielectronic']
            abundance = chdata.Abundance[self.AbundanceName]['abundance'][Z - 1]
            if verbose:
                print(' %5i %5s abundance = %10.2e '%(Z, const.El[Z-1],  abundance))
            if verbose:
                print(' doing ion %s for the following processes %s'%(akey, self.Todo[akey]))
            if 'ff' in self.Todo[akey]:
                if verbose:
                    print(' calculating ff continuum for :  %s'%(akey))
                FF = chianti.core.continuum(akey, temperature, abundanceName=self.AbundanceName, em=em)
                FF.freeFree(wavelength)
#                if nTempDen == 1:
#                    freeFree += cont.FreeFree['rate']*em[0]
#                else:
#                    for iTempDen in range(nTempDen):
#                        freeFree[iTempDen] += cont.FreeFree['rate'][iTempDen]*em[iTempDen]
                freeFree += FF.FreeFree['rate']
            if 'fb' in self.Todo[akey]:
                if verbose:
                    print(' calculating fb continuum for :  %s'%(akey))
                try:
                    FB = chianti.core.continuum(akey, temperature, abundanceName=self.AbundanceName, em=em)
                    FB.freeBound(wavelength)
                    if 'errorMessage' not in list(cont.FreeBound.keys()):
                        #  an fblvl file exists for this ions
                        freeBound += FB.FreeBound['rate']
                except:
                    cont = chianti.core.continuum(akey, temperature, abundanceName=self.AbundanceName)
                    cont.freeBound(wavelength)
#                if 'errorMessage' not in list(cont.FreeBound.keys()):
#                    #  an fblvl file exists for this ions
#                    if nTempDen == 1:
#                        freeBound += cont.FreeBound['rate']*em[0]
#                    else:
#                        for iTempDen in range(nTempDen):
#                            freeBound[iTempDen] += cont.FreeBound['rate'][iTempDen]*em[iTempDen]
#                        freeBound += cont.FreeBound['rate']
            if 'line' in self.Todo[akey]:
                if verbose:
                    print(' calculating spectrum for  :  %s'%(akey))
                thisIon = chianti.core.ion(akey, temperature, eDensity, abundanceName=self.AbundanceName)
                thisIon.intensity(wvlRange=wvlRange, allLines=allLines, em=em)
                self.IonsCalculated.append(akey)
                if 'errorMessage' not in  list(thisIon.Intensity.keys()):
                    self.Finished.append(akey)
                    thisIon.spectrum(wavelength, filter=filter, allLines=allLines, em=em)
                    if keepIons:
                        self.IonInstances[akey] = copy.deepcopy(thisIon)
                    if setupIntensity:
                        for bkey in self.Intensity:
                            self.Intensity[bkey] = np.hstack((copy.copy(self.Intensity[bkey]), thisIon.Intensity[bkey]))
                    else:
                        setupIntensity = 1
                        self.Intensity  = thisIon.Intensity
                    lineSpectrum += thisIon.Spectrum['intensity']
#                            if nTempDen == 1:
#                                lineSpectrum += thisIon.Spectrum['intensity']
#                            else:
#                                for iTempDen in range(nTempDen):
#                                    lineSpectrum[iTempDen] += thisIon.Spectrum['intensity'][iTempDen]
                else:
                    if verbose:
                        print(thisIon.Intensity['errorMessage'])
                # get 2 photon emission for H and He sequences
                if (Z - ionstage) in [0, 1] and not dielectronic:
                    thisIon.twoPhoton(wavelength)
                    twoPhoton += thisIon.TwoPhoton['rate']

        self.FreeFree = {'wavelength':wavelength, 'intensity':freeFree.squeeze()}
        self.FreeBound = {'wavelength':wavelength, 'intensity':freeBound.squeeze()}
        self.LineSpectrum = {'wavelength':wavelength, 'intensity':lineSpectrum.squeeze()}
        self.TwoPhoton = {'wavelength':wavelength, 'intensity':twoPhoton.squeeze()}
        #
        #
        total = freeFree + freeBound + lineSpectrum + twoPhoton
        self.Total = total
        t2 = datetime.now()
        dt=t2-t1
        print(' elapsed seconds = %12.3f'%(dt.seconds))
        if nTempDen == 1:
            integrated = total
        else:
            integrated = total.sum(axis=0)
        #
        if type(label) == type(''):
            if hasattr(self, 'Spectrum'):
                self.Spectrum[label] = {'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':em, 'ions':self.IonsCalculated,
                'Abundance':self.AbundanceName, 'xlabel':xlabel, 'ylabel':ylabel}
            else:
                self.Spectrum = {label:{'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':em, 'ions':self.IonsCalculated,
                'Abundance':self.AbundanceName, 'xlabel':xlabel, 'ylabel':ylabel}}
        else:
            self.Spectrum ={'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated,  'ions':self.IonsCalculated,
            'Abundance':self.AbundanceName, 'xlabel':xlabel, 'ylabel':ylabel}
        #
        # -----------------------------------------------------------------------
        #
#        self.Spectrum ={'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'ions':self.IonsCalculated, 'Abundance':self.AbundanceName}
    #
    # ----------------------------------------------------------------------------------------------
    #
class bunch(_ionTrails, _specTrails):
    '''
    Calculate the emission line spectrum as a function of temperature and density.

    'bunch' is very similar to 'spectrum' except that continuum is not calculated and the spectrum
    is not convolved over a filter.  However, this can be done with the inherited convolve method

    one of the convenient things is that all of the instantiated ion classes, determined through such keywords as 'elementList',
    'ionList', and 'minAbund' are kept in a dictionary self.IonInstances where self.IonInstances['mg_7'] is the class instance of
    chianti.core.ion for 'mg_7'.  All its methods and attributes are available.

    includes elemental abundances and ionization equilibria

    the set of abundances, a file in $XUVTOP/abundance, can be set with the keyword argument 'abundanceName'
    temperature and density can be arrays but, unless the size of either is one (1),
    the two must have the same size

    Inherited methods include 'intensityList', 'intensityRatio' (between lines of different ions), and 'intensityRatioSave'
    and 'convolve'.

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
    temperature/density
    '''
    #
    # ------------------------------------------------------------------------------------
    #
    def __init__(self, temperature, eDensity, wvlRange, elementList=0, ionList=0, minAbund=0, keepIons=0, em=0, abundanceName=0, verbose=0, allLines=1):
        #
        t1 = datetime.now()
        # creates Intensity dict from first ion calculated
        setupIntensity = 0
        #
        self.Defaults=defaults
        temperature = np.asarray(temperature, 'float64')
        self.Temperature = temperature
        eDensity = np.asarray(eDensity, 'float64')
        self.EDensity = eDensity

        #
        self.EDensity = np.asarray(eDensity,'float64')
        self.NEDens = self.EDensity.size
        ndens = self.EDensity.size
        ntemp = self.Temperature.size
        tst1 = ndens == ntemp
        tst1a = ndens != ntemp
        tst2 = ntemp > 1
        tst3 = ndens > 1
        tst4 = ndens > 1 and ntemp > 1
        if tst1 and ntemp == 1:
            self.NTempDen = 1
        elif tst1a and (tst2 or tst3) and not tst4:
            self.NTempDen = ntemp*ndens
            if ntemp == self.NTempDen and ndens != self.NTempDen:
                self.EDensity = np.ones_like(self.Temperature)*self.EDensity
            elif ndens == self.NTempDen and ntemp != self.NTempDen:
                self.Temperature = np.ones_like(self.EDensity)*self.Temperature
        elif tst1 and tst4:
            self.NTempDen = ntemp
        #
        if type(em) == int and em == 0:
            em = np.ones(self.NTempDen, 'float64')
        elif type(em) == float and em > 0.:
            em = np.ones(self.NTempDen, 'float64')*em
        elif type(em) == list or type(em) == tuple:
            em = np.asarray(em, 'float64')
        self.Em = em
        #
        #if em != 0:
            #if type(em) == type(float):
                #if nTempDen > 1:
                    #em = np.ones_like(self.Temperature)*em
                    #nEm = nTempDen
                #else:
                    #nEm = 1
            #else:
                #em = np.asarray(em, 'float64')
                #nEm = em.size
                #if nEm != nTempDen:
                    #print(' the emission measure array must be the same size as the temperature/density array')
                    #return
        #else:
            #nEm = 0
        #self.NEm = nEm
        #self.AbundanceName = defaults['abundfile']
        #self.AbundanceAll = chdata.AbundanceAll
        #
        if abundanceName:
            if abundanceName in list(chdata.Abundance.keys()):
                self.AbundanceName = abundanceName
            else:
                abundChoices = list(chdata.Abundance.keys())
#                for one in wvl[topLines]:
#                    wvlChoices.append('%12.3f'%(one))
                abundChoice = chGui.gui.selectorDialog(abundChoices,label='Select Abundance name')
                abundChoice_idx = abundChoice.selectedIndex
                self.AbundanceName = abundChoices[abundChoice_idx[0]]
                abundanceName = self.AbundanceName
                print((' Abundance chosen:  %s '%(self.AbundanceName)))
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
        ionInfo = chio.masterListInfo()
        #        #
#        self.Intensity = {'ionS':[], 'lvl1':[], 'lvl2':[], 'wvl':np.ndarray, 'pretty1':np.ndarray, 'pretty2':np.ndarray, 'intensity':np.zeros((nTempDen, 0),'float64'), 'obs':np.ndarray }
        self.IonsCalculated = []
        if keepIons:
            self.IonInstances = {}
        self.Finished = []
        #
        # also needed by ionGate
        self.WvlRange = np.asarray(wvlRange, 'float64')
        #
        self.ionGate(elementList = elementList, ionList = ionList, minAbund=minAbund, doContinuum=0, verbose = verbose)
        #
        for ionS in sorted(self.Todo.keys()):
            nameStuff = util.convertName(ionS)
            Z = nameStuff['Z']

            if verbose:
                print(' calculating %s'%(ionS))
            thisIon = chianti.core.ion(ionS, temperature, eDensity, abundance=abundAll[Z-1])
            thisIon.intensity(wvlRange=wvlRange, allLines = allLines,  em=em)
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
        return
    #
    # -------------------------------------------------------------------------
    #
