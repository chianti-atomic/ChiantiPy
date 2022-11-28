from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
np.seterr(over='ignore')

from .Continuum import continuum
from .Ion import ion
import ChiantiPy.tools.data as chdata
import ChiantiPy.tools.constants as const
import ChiantiPy.tools.util as util
import ChiantiPy.tools.io as chio
from ChiantiPy.base import specTrails
from ChiantiPy.base import ionTrails


class radLoss(ionTrails, specTrails):
    '''
    Calculate the radiative emission loss rate as a function of temperature and density.

    includes elemental abundances or ionization equilibria

    temperature and density can be arrays but, unless the size of either is one (1),
    the two must have the same size


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

    abundance: to select a particular set of abundances, set abundance to the name of a CHIANTI abundance file,
        without the '.abund' suffix, e.g. 'sun_photospheric_1998_grevesse'
        If set to a blank (''), a gui selection menu will popup and allow the
        selection of an set of abundances
    '''
    def __init__(self, temperature, eDensity, elementList=None, ionList = None, minAbund=None,
        doContinuum=True, doLines=True, abundance=None, verbose=0, allLines=1):
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
        self.Defaults=chdata.Defaults
        self.argCheck(temperature=temperature, eDensity=eDensity, pDensity=None, em=None)
        self.Labels = util.units(chdata.Defaults)

        #
        if abundance is not None:
            ab = chio.abundanceRead(abundance)
            abundAll = ab['abundance']
            self.AbundanceName = abundance
        else:
            self.AbundanceName = self.Defaults['abundfile']
        #
            abundAll = chdata.Abundance[self.AbundanceName]['abundance']
        #
#        abundAll = chdata.Abundance[self.AbundanceName]['abundance']
#        # needed by ionGate
        self.AbundAll = abundAll
        self.Abundance = abundAll
        #
#        nonzed = self.Abundance > 0.
#        minAbundAll = self.Abundance[nonzed].min()
        # if minAbund is even set
#        if minAbund:
#            if minAbund < minAbundAll:
#                minAbund = minAbundAll
#        self.MinAbund = minAbund
#        ionInfo = util.masterListInfo()
        #

        freeFreeLoss = np.zeros_like(self.Temperature)
        freeBoundLoss = np.zeros_like(self.Temperature)
        twoPhotonLoss = np.zeros_like(self.Temperature)
        boundBoundLoss = np.zeros_like(self.Temperature)
        twoPhotonLoss = np.zeros_like(self.Temperature)
        #
        self.IonsCalculated = []

        self.Finished = []
        #
        self.WvlRange = [0., 1.e+30]
        #
        self.ionGate(elementList = elementList, ionList = ionList, minAbund=minAbund, doLines=doLines,
            doContinuum=doContinuum, doWvlTest=0, verbose=False)
        #
        #
        for akey in sorted(self.Todo.keys()):
            zStuff = util.convertName(akey)
            Z = zStuff['Z']
            ionstage = zStuff['Ion']
            dielectronic = zStuff['Dielectronic']
            abundanceZ = self.Abundance[Z - 1]
            if verbose:
                print(' %5i %5s abundance = %10.2e '%(Z, const.El[Z-1],  abundanceZ))
            if verbose:
                print(' doing ion %s for the following processes %s'%(akey, self.Todo[akey]))
            if ionstage != 1:
                if verbose:
                    print(' calculating ff continuum for :  %s'%(akey))
                if 'ff' in self.Todo[akey]:
                    # need to skip the neutral
                    cont = continuum(akey, temperature, abundance=abundanceZ)
                    cont.freeFreeLoss()
                    freeFreeLoss += cont.FreeFreeLoss['rate']
                if 'fb' in self.Todo[akey]:
                    if verbose:
                        print(' calculating fb continuum for :  %s'%(akey))
                    cont = continuum(akey, temperature, abundance=abundanceZ)
                    cont.freeBoundLoss(verbose=verbose)
                    if 'errorMessage' not in cont.FreeBoundLoss.keys():
                        freeBoundLoss += cont.FreeBoundLoss['rate']
            if 'line' in self.Todo[akey]:
                if verbose:
                    print(' calculating spectrum for  :  %s'%(akey))
                thisIon = ion(akey, temperature, eDensity, abundance=abundanceZ)
                thisIon.intensity(allLines=allLines)
                self.IonsCalculated.append(akey)
                if 'errorMessage' not in  thisIon.Intensity.keys():
                    self.Finished.append(akey)
                    thisIon.boundBoundLoss()
                    boundBoundLoss += thisIon.BoundBoundLoss['rate']
                else:
                    if verbose:
                        print(thisIon.Intensity['errorMessage'])
                # get 2 photon emission for H and He sequences
                if (Z - ionstage) in [0, 1] and not dielectronic:
                    thisIon.twoPhotonLoss()
                    twoPhotonLoss += thisIon.TwoPhotonLoss['rate']
        self.FreeFreeLoss = freeFreeLoss
        self.FreeBoundLoss = freeBoundLoss
        self.BoundBoundLoss = boundBoundLoss
        self.TwoPhotonLoss = twoPhotonLoss
        #
        total = freeFreeLoss + freeBoundLoss + boundBoundLoss + twoPhotonLoss
        t2 = datetime.now()
        dt=t2-t1
        print(' elapsed seconds = %10.2e'%(dt.seconds))
        xlabel = self.Labels['radlossTlabel']
        ylabel = self.Labels['radlossYlabel']
        self.RadLoss = {'rate':total, 'temperature':self.Temperature, 'density':self.EDensity,
            'minAbund':minAbund, 'abundance':self.AbundanceName, 'ylabel':ylabel, 'xlabel':xlabel}
    #
    # -------------------------------------------------------------------
    #
    def radLossPlot(self, doTitle=False):
        '''
        to plot the radiative losses vs temperature

        Parameters
        ----------

        doTitle:  `bool`

            if True, a title is applied to the plot.  The default is for no title

        '''
        fontsize = 16
        temp = self.RadLoss['temperature']
        rate = self.RadLoss['rate']
        plt.loglog(temp, rate,  'k',  lw=2)
        plt.xlabel(self.RadLoss['xlabel'],fontsize=fontsize)
        plt.ylabel(self.RadLoss['ylabel'],fontsize=fontsize)
        plt.xlim(left = temp.min(),  right = temp.max())
        if doTitle:
            if hasattr(self, 'AbundanceName'):
                title = 'Radiative loss rate,  %s'%(self.AbundanceName)
            if hasattr(self, 'MinAbund'):
                title += ' minAbund = %10.2e'%(self.MinAbund)
            if self.EDensity.size == 1:
                title += ', density = %10.2e'%(self.EDensity)
            plt.title(title, fontsize=fontsize)
        plt.tight_layout()
