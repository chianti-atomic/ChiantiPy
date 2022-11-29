from datetime import datetime
import copy
import multiprocessing as mp

import numpy as np
import matplotlib.pyplot as plt
np.seterr(over='ignore')

#from .Continuum import continuum
#from .Ion import ion
import ChiantiPy.tools.data as chdata
import ChiantiPy.tools.constants as const
import ChiantiPy.tools.util as util
import ChiantiPy.tools.io as chio
from ChiantiPy.base import specTrails
from ChiantiPy.base import ionTrails
import ChiantiPy.tools.mputil as mputil


class mradLoss(ionTrails, specTrails):
    '''
    Calculate the radiative emission loss rate as a function of temperature and density.

    this is the multiprocessing version of radloss

    includes elemental abundances or ionization equilibria

    temperature and density can be arrays but, unless the size of either is one (1),
    the two must have the same size


    A selection of ions can be make with ionList containing the names of
    the desired lines in Chianti notation,
    i.e. C VI = c_6

    a minimum abundance can be specified so that the calculation can be speeded up by excluding
    elements with a low abundance. With solar photospheric abundances

    setting minAbund = 1.e-4 will include H, He, C, O, Ne
    setting minAbund = 2.e-5 adds  N, Mg, Si, S, Fe
    setting minAbund = 1.e-6 adds  Na, Al, Ar, Ca, Ni

    Setting em will multiply the spectrum at each temperature by the value of em.

    em [for emission measure], can be a float or an array of the same length as the
    temperature/density.

    abundance: to select a particular set of abundances, set abundance to the name of a
        CHIANTI abundance file, without the '.abund' suffix, e.g. 'sun_photospheric_1998_grevesse'
        If set to a blank (''), a gui selection menu will popup and allow the
        selection of an set of abundances

    Parameters
    --------------

    temperature: `float`, `list`, `ndarray`
        the temperature(s) in K

    eDensity: float, ndarray
        eDensity: electron density in :math:`\mathrm{cm^{-3}}`

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

    abundance:  `str`
        abuncance:  the file name of the abuncance set to be used
            must be one in the $XUVTOP/abund directory

    allLInes:  `bool`
        allLines:  whether or not to include unobserved lines

    verbose:  `bool`
        verbose:  whether to allow certain print statements

    '''
    def __init__(self, temperature, eDensity, elementList=None, ionList = None, minAbund=None,
        doContinuum=True, doLines=True, abundance=None, verbose=0, allLines=1,  proc=4):

        timeout = 0.01
        t1 = datetime.now()
        masterlist = chdata.MasterList
        # use the ionList but make sure the ions are in the database

        setupIntensity = False
        #
        self.Defaults = chdata.Defaults
        self.Labels = util.units(chdata.Defaults)

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

        self.AllLines = allLines
        self.MinAbund = minAbund
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
#        nTempDens = self.NTempDens
        freeFreeLoss = np.zeros_like(self.Temperature)
        freeBoundLoss = np.zeros_like(self.Temperature)
        twoPhotonLoss = np.zeros_like(self.Temperature)
        boundBoundLoss = np.zeros_like(self.Temperature)
        twoPhotonLoss = np.zeros_like(self.Temperature)
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
        self.Finished = []
        #
#        self.WvlRange = [0., 1.e+30]
        #
        self.ionGate(elementList = elementList, ionList = ionList, minAbund=minAbund,
            doContinuum=doContinuum, doLines=doLines, doWvlTest=0, verbose=False)
        #
        for akey in sorted(self.Todo.keys()):
            zStuff = util.convertName(akey)
            Z = zStuff['Z']
            abundance = self.Abundance[Z - 1]
#            if verbose:
#                print(' %5i %5s abundance = %10.2e '%(Z, const.El[Z-1],  abundance))
            if verbose:
                if self.Todo[akey] != '':
                    print(' doing ion %s for the following processes %s'%(akey, self.Todo[akey]))
            if 'ff' in self.Todo[akey]:
                ffWorkerQ.put((akey, temperature, abundance))
            if 'fb' in self.Todo[akey]:
                fbWorkerQ.put((akey, temperature, abundance))
            if 'line' in self.Todo[akey]:
                ionWorkerQ.put((akey, temperature, eDensity, allLines, abundance, doContinuum))
        #
        ffWorkerQSize = ffWorkerQ.qsize()
        fbWorkerQSize = fbWorkerQ.qsize()
        ionWorkerQSize = ionWorkerQ.qsize()

        if doContinuum:
            ffProcesses = []
            for i in range(proc):
                p = mp.Process(target=mputil.doFfLossQ, args=(ffWorkerQ, ffDoneQ))
                p.start()
                ffProcesses.append(p)
    #       timeout is not necessary
            for p in ffProcesses:
                if p.is_alive():
                    p.join(timeout=timeout)
            #
            for iff in range(ffWorkerQSize):
                thisFreeFree = ffDoneQ.get()
                freeFreeLoss += thisFreeFree['rate']
            for p in ffProcesses:
                if not isinstance(p, str):
                    p.terminate()
        #
            fbProcesses = []
            for i in range(proc):
                p = mp.Process(target=mputil.doFbLossQ, args=(fbWorkerQ, fbDoneQ))
                p.start()
                fbProcesses.append(p)
    #       timeout is not necessary
            for p in fbProcesses:
                if p.is_alive():
                    p.join(timeout=timeout)
            #
            for ifb in range(fbWorkerQSize):
                thisFreeBound = fbDoneQ.get()
                if 'errorMessage' not in thisFreeBound.keys():
                    freeBoundLoss += thisFreeBound['rate'].squeeze()

            for p in fbProcesses:
                if not isinstance(p, str):
                    p.terminate()
        #
        if doLines:
            ionProcesses = []
            if ionWorkerQSize < proc:
                proc = ionWorkerQSize
            for i in range(proc):
                p = mp.Process(target=mputil.doIonLossQ, args=(ionWorkerQ, ionDoneQ))
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
                thisRadLoss = thisIon.BoundBoundLoss
                if not 'errorMessage' in sorted(thisRadLoss.keys()):
                    self.Finished.append(ionS)
                    boundBoundLoss += thisRadLoss['rate']
#                    if setupIntensity:
#                        for akey in sorted(self.Intensity.keys()):
#                            self.Intensity[akey] = np.hstack((copy.copy(self.Intensity[akey]),
#                                thisIntensity[akey]))
#                    else:
#                        setupIntensity = True
#                        self.Intensity  = thisIntensity
                    #
#                    if not 'errorMessage' in sorted(thisIon.Spectrum.keys()):
                    boundBoundLoss += thisIon.BoundBoundLoss['rate']
                   # check for two-photon emission
                    if len(out) == 3:
                        tp = out[2]
                        twoPhotonLoss += tp['rate'].squeeze()
                else:
                    if 'errorMessage' in sorted(thisRadLoss.keys()):
                        print(thisRadLoss['errorMessage'])
                #
            for p in ionProcesses:
                if not isinstance(p, str):
                    p.terminate()

#        for akey in sorted(self.Todo.keys()):
#            zStuff = util.convertName(akey)
#            Z = zStuff['Z']
#            ionstage = zStuff['Ion']
#            dielectronic = zStuff['Dielectronic']
#            abundanceZ = self.Abundance[Z - 1]
#            if verbose:
#                print(' %5i %5s abundance = %10.2e '%(Z, const.El[Z-1],  abundanceZ))
#            if verbose:
#                print(' doing ion %s for the following processes %s'%(akey, self.Todo[akey]))
#            if ionstage != 1:
#                if verbose:
#                    print(' calculating ff continuum for :  %s'%(akey))
#                if 'ff' in self.Todo[akey]:
#                    # need to skip the neutral
#                    cont = continuum(akey, temperature, abundance=abundanceZ)
#                    cont.freeFreeLoss()
#                    freeFreeLoss += cont.FreeFreeLoss['rate']
#                if 'fb' in self.Todo[akey]:
#                    if verbose:
#                        print(' calculating fb continuum for :  %s'%(akey))
#                    cont = continuum(akey, temperature, abundance=abundanceZ)
#                    cont.freeBoundLoss(verbose=verbose)
#                    if 'errorMessage' not in cont.FreeBoundLoss.keys():
#                        freeBoundLoss += cont.FreeBoundLoss['rate']
#            if 'line' in self.Todo[akey]:
#                if verbose:
#                    print(' calculating spectrum for  :  %s'%(akey))
#                thisIon = ion(akey, temperature, eDensity, abundance=abundanceZ)
#                thisIon.intensity(allLines=allLines)
#                self.IonsCalculated.append(akey)
#                if 'errorMessage' not in  thisIon.Intensity.keys():
#                    self.Finished.append(akey)
#                    thisIon.boundBoundLoss()
#                    boundBoundLoss += thisIon.BoundBoundLoss['rate']
#                else:
#                    if verbose:
#                        print(thisIon.Intensity['errorMessage'])
#                # get 2 photon emission for H and He sequences
#                if (Z - ionstage) in [0, 1] and not dielectronic:
#                    thisIon.twoPhotonLoss()
#                    twoPhotonLoss += thisIon.TwoPhotonLoss['rate']
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
    def radLossPlot(self, title=0):
        '''
        to plot the radiative losses vs temperature
        '''
        fontsize = 16
        temp = self.RadLoss['temperature']
        rate = self.RadLoss['rate']
        plt.loglog(temp, rate)
#        plt.ylabel(r'erg  s$^{-1}$  ($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$',fontsize=fontsize)
        plt.xlabel(self.RadLoss['xlabel'],fontsize=fontsize)
        plt.ylabel(self.RadLoss['ylabel'],fontsize=fontsize)
        if title:
            title = 'Radiative loss rate,  minAbund = %10.2e'%(self.MinAbund)
            if self.EDensity.size == 1:
                title += ', density = %10.2e'%(self.EDensity)
            plt.title(title, fontsize=fontsize)
        plt.tight_layout()
