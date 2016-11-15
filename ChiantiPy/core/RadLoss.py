from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
np.seterr(over='ignore')
#import ChiantiPy.core as ch
import ChiantiPy
import ChiantiPy.tools.data as chdata
import ChiantiPy.tools.constants as const
import ChiantiPy.tools.util as util
import ChiantiPy.Gui as chgui
from ._SpecTrails import specTrails
#
defaults = chdata.Defaults
#chInteractive = chdata.chInteractive
#if chInteractive:
#    import pylab as pl
#else:
#    import matplotlib
#    matplotlib.use('Agg')
#    import matplotlib.pyplot as pl

class radLoss(specTrails):
    '''
    Calculate the emission spectrum as a function of temperature and density.

    includes elemental abundances or ionization equilibria

    temperature and density can be arrays but, unless the size of either is one (1),
    the two must have the same size

    the returned spectrum will be convolved with a filter of the specified width on the
    specified wavelength array

    the default filter is gaussianR with a resolving power of 1000.  Other filters,
    such as gaussian, box and lorentz, are available in ChiantiPy.filters.  When using the box filter,
    the width should equal the wavelength interval to keep the units of the continuum and line
    spectrum the same.

    A selection of ions can be make with ionList containing the names of
    the desired lines in Chianti notation, i.e. C VI = c_6

    a minimum abundance can be specified so that the calculation can be speeded up by excluding
    elements with a low abundance. With solar photospheric abundances -

    setting minAbund = 1.e-4 will include H, He, C, O, Ne
    setting minAbund = 2.e-5 adds  N, Mg, Si, S, Fe
    setting minAbund = 1.e-6 adds  Na, Al, Ar, Ca, Ni

    Setting em will multiply the spectrum at each temperature by the value of em.

    em [for emission measure], can be a float or an array of the same length as the
    temperature/density.
    '''
    def __init__(self, temperature, eDensity, elementList=0, ionList = 0, minAbund=0, doContinuum=1, abundanceName=0, verbose=0, allLines=1, keepIons=0):
        t1 = datetime.now()
        masterlist = chdata.MasterList
        # use the ionList but make sure the ions are in the database
        if ionList:
            alist=[]
            for one in ionList:
                if masterlist.count(one):
                    alist.append(one)
                else:
                    if verbose:
                        pstring = ' %s not in CHIANTI database'%(one)
                        print(pstring)
            masterlist = alist
        self.Defaults=defaults
        self.Temperature = np.asarray(temperature, 'float64')
        nTemp = self.Temperature.size
        self.EDensity = np.asarray(eDensity, 'float64')
        nDen = self.EDensity.size
        nTempDen = max([nTemp, nDen])

        if not abundanceName:
            self.AbundanceName = self.Defaults['abundfile']
        else:
            if abundanceName in sorted(chdata.Abundance.keys()):
                self.AbundanceName = abundanceName
            else:
                abundChoices = sorted(chdata.Abundance.keys())
#                for one in wvl[topLines]:
#                    wvlChoices.append('%12.3f'%(one))
                abundChoice = chgui.gui.selectorDialog(abundChoices,label='Select Abundance name')
                abundChoice_idx = abundChoice.selectedIndex
                self.AbundanceName = abundChoices[abundChoice_idx[0]]
                abund = self.AbundanceName
                print(' Abundance chosen:  %s '%(self.AbundanceName))
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
#        ionInfo = util.masterListInfo()
        #
        freeFreeLoss = np.zeros((nTempDen), 'float64').squeeze()
        freeBoundLoss = np.zeros((nTempDen), 'float64').squeeze()
        twoPhotonLoss = np.zeros((nTempDen), 'float64').squeeze()
        boundBoundLoss = np.zeros((nTempDen), 'float64').squeeze()
        twoPhotonLoss = np.zeros((nTempDen), 'float64').squeeze()
        #
        self.IonsCalculated = []
        if keepIons:
            self.IonInstances = {}
        self.Finished = []
        #
        self.WvlRange = [0., 1.e+30]
        #
        self.ionGate(elementList = elementList, ionList = ionList, minAbund=minAbund, doContinuum=doContinuum, doWvlTest=0, verbose=verbose)
        #
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
            if ionstage != 1:
                if verbose:
                    print(' calculating ff continuum for :  %s'%(akey))
                if 'ff' in self.Todo[akey]:
                    # need to skip the neutral
                        cont = ChiantiPy.core.continuum(akey, temperature, abundanceName=self.AbundanceName)
                        cont.freeFreeLoss()
                        freeFreeLoss += cont.FreeFreeLoss['rate']
    #                if nTempDen == 1:
    #                    freeFree += cont.FreeFreeLoss['rate']
    #                else:
    #                    for iTempDen in range(nTempDen):
    #                        freeFree[iTempDen] += cont.FreeFree['rate'][iTempDen]*em[iTempDen]
    #                freeFree += cont.FreeFree['rate']
                if 'fb' in self.Todo[akey]:
                    if verbose:
                        print(' calculating fb continuum for :  %s'%(akey))
    #                try:
                    # does cont already exist - i.e did we create it for ff
                    if hasattr(cont, 'FreeFreeLoss'):
                        cont.freeBoundLoss()
                    else:
                        cont = ChiantiPy.core.continuum(akey, temperature, abundanceName=self.AbundanceName)
                        cont.freeBoundLoss()
                    if 'errorMessage' not in list(cont.FreeBoundLoss.keys()):
                        #  an fblvl file exists for this ions
                        freeBoundLoss += cont.FreeBoundLoss['rate']
#                except:
#                    cont = ChiantiPy.core.continuum(akey, temperature, abundanceName=self.AbundanceName)
#                    cont.freeBoundLoss()
#                if 'errorMessage' not in list(cont.FreeBound.keys()):
#                    #  an fblvl file exists for this ions
#                    freeBoundLoss += cont.FreeBoundLoss['rate']
#                    if nTempDen == 1:
#                        freeBoundLoss += cont.FreeBoundLoss['rate']
#                    else:
#                        for iTempDen in range(nTempDen):
#                            freeBound[iTempDen] += cont.FreeBound['rate'][iTempDen]*em[iTempDen]
#                        freeBound += cont.FreeBound['rate']
            if 'line' in self.Todo[akey]:
                if verbose:
                    print(' calculating spectrum for  :  %s'%(akey))
                thisIon = ChiantiPy.core.ion(akey, temperature, eDensity, abundanceName=self.AbundanceName)
                thisIon.intensity(allLines=allLines)
                self.IonsCalculated.append(akey)
                if 'errorMessage' not in  list(thisIon.Intensity.keys()):
                    self.Finished.append(akey)
                    thisIon.boundBoundLoss()
#                    self.IonInstances[akey] = copy.deepcopy(thisIon)
#                    if setupIntensity:
#                        for bkey in self.Intensity:
#                            self.Intensity[bkey] = np.hstack((copy.copy(self.Intensity[bkey]), thisIon.Intensity[bkey]))
#                    else:
#                        setupIntensity = 1
#                        self.Intensity  = thisIon.Intensity
                    boundBoundLoss += thisIon.BoundBoundLoss['rate']
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
                    thisIon.twoPhotonLoss()
                    twoPhotonLoss += thisIon.TwoPhotonLoss['rate']

#        for iz in range(31):
#            abundance = chdata.Abundance['abundance'][iz-1]
#            if abundance >= minAbund:
#                print(' %5i %5s abundance = %10.2e '%(iz, const.El[iz-1],  abundance))
#                #
#                for ionstage in range(1, iz+2):
#                    ionS = util.zion2name(iz, ionstage)
##                   print ' ionS = ', ionS
#                    masterListTest = ionS in masterlist
#                    masterListInfoTest = ionS in sorted(ionInfo.keys())
#                    if masterListTest or masterListInfoTest:
#                        ioneqTest = (self.Temperature.max() >= ionInfo[ionS]['tmin']) and (self.Temperature.min() <= ionInfo[ionS]['tmax'])
#                    # construct similar test for the dielectronic files
#                    ionSd = util.zion2name(iz, ionstage, dielectronic=1)
#                    masterListTestD = ionSd in masterlist
#                    masterListInfoTestD = ionSd in sorted(ionInfo.keys())
#                    if masterListTestD or masterListInfoTestD:
#                        ioneqTestD = (self.Temperature.max() >= ionInfo[ionSd]['tmin']) and (self.Temperature.min() <=ionInfo[ionSd]['tmax'])
#                    ionstageTest = ionstage > 1
#                    if ionstageTest and ioneqTest and doContinuum:
#                        # ionS is the target ion, cannot be the neutral for the continuum
#                        print(' calculating continuum for %s'%(ionS))
#                        cont = ChiantiPy.core.continuum(ionS, temperature, abund=abund)
#                        cont.freeFreeLoss()
#    #                   print dir(thisIon)
#    #                   print ' wvl = ', thisIon.FreeFree['wvl']
##                        if nTempDen ==1:
#                        freeFreeLoss += cont.FreeFreeLoss['rate']
##                        else:
##                            for iTempDen in range(nTempDen):
##                                freeFreeLoss[iTempDen] += cont.FreeFreeLoss['rate'][iTempDen]
#                    #
#                        cont.freeBoundLoss()
#                        if 'errorMessage' not in sorted(cont.FreeBoundLoss.keys()):
#                            #  an fblvl file exists for this ions
##                            if nTempDen == 1:
#                            freeBoundLoss += cont.FreeBoundLoss['rate']
##                            else:
##                                freeBound[iTempDen] += cont.FreeBound['rate'][iTempDen]
#                    if masterListTest and ioneqTest:
#                        print(' calculating spectrum for  %s  '%(ionS))
#                        thisIon = ChiantiPy.core.ion(ionS, temperature, density, abund=abund)
##                       print ' dir = ', dir(thisIon)
##                        thisIon.emiss(wvlRange = wvlRange, allLines=allLines)
#                        thisIon.boundBoundLoss( allLines=allLines)
#                        # check that there are lines in this wavelength range
#                        if 'errorMessage' not in  sorted(thisIon.BoundBoundLoss.keys()):
#                            thisIon.boundBoundLoss()
##                           intensity = thisIon.Intensity['intensity']
##                            if nTempDen == 1:
#                            boundBoundLoss += thisIon.BoundBoundLoss['rate']
##                            else:
##                                for iTempDen in range(nTempDen):
##                                    lineSpectrum[iTempDen] += thisIon.Spectrum['intensity'][iTempDen]
#                        # get 2 photon emission for H and He sequences
#                        if (iz - ionstage) in [0, 1]:
#                            thisIon.twoPhotonLoss()
#                            twoPhotonLoss += thisIon.TwoPhotonLoss['rate']
#                    # get dielectronic lines
#                    if masterListTestD and ioneqTestD:
#                        print(' calculating spectrum for  %s '%(ionSd))
#                        thisIon = ChiantiPy.core.ion(ionSd, temperature, density, abund=abund)
##                       print ' dir = ', dir(thisIon)
##                       have to do all lines for the dielectronic satellites
##                        thisIon.emiss(allLines=1)
#                        thisIon.intensity(allLines=allLines)
#                        # check that there are lines in this wavelength range - probably not redundant
#                        if 'errorMessage' not in  sorted(thisIon.Intensity.keys()):
#                            thisIon.boundBoundLoss()
##                            if nTempDen == 1:
#                            boundBoundLoss += thisIon.BoundBoundLoss['rate']
##                            else:
##                                for iTempDen in range(nTempDen):
##                                    lineSpectrum[iTempDen] += thisIon.Spectrum['intensity'][iTempDen]
        self.FreeFreeLoss = freeFreeLoss
        self.FreeBoundLoss = freeBoundLoss
        self.BoundBoundLoss = boundBoundLoss
        self.TwoPhotonLoss = twoPhotonLoss
        #
        total = freeFreeLoss + freeBoundLoss + boundBoundLoss + twoPhotonLoss
        t2 = datetime.now()
        dt=t2-t1
        print(' elapsed seconds = %10.2e'%(dt.seconds))
        self.RadLoss = {'rate':total, 'temperature':self.Temperature, 'density':self.EDensity, 'minAbund':minAbund, 'abundance':self.AbundanceName}
    #
    # -------------------------------------------------------------------
    #
    def radLossPlot(self):
        '''
        to plot the radiative losses vs temperature
        '''
        fontsize = 16
        temp = self.RadLoss['temperature']
        rate = self.RadLoss['rate']
        pl.loglog(temp, rate)
        plt.xlabel('Temperature (K)',fontsize=fontsize)
        plt.ylabel(r'erg  s$^{-1}$  ($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$',fontsize=fontsize)
#            title = 'Radiative loss rate,  minAbund = %6.1e'%(self.MinAbund)
#            if self.Density.size == 1:
#                title += ', density = %6.1e'%(self.Density)
#            plt.title(title, fontsize=fontsize)
