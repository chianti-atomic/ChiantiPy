'''
classes and methods to analyze observations of astrophysical spectra
'''
import os
from datetime import datetime
import pickle

try:
    import multiprocessing as mp
except:
    print(' your version of Python does not support multiprocessing \n you will not be able to use mgofnt')
#
import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import ChiantiPy.core as ch
import ChiantiPy.tools.io as io
import ChiantiPy.tools.util as util
import ChiantiPy.tools.constants as const
import ChiantiPy.tools.data as chdata
import ChiantiPy.Gui as chGui
from ChiantiPy.base import ionTrails    #
    # --------------------------------------------------------------------------
    #
def doDemGofntQ(inQueue, outQueue):
    '''
    helper for multiprocessing with maker.mgofnt()
    '''
    for inputs in iter(inQueue.get, 'STOP'):
        ionS = inputs[0]
#        print(' helper doing '+ionS)
        temperature = inputs[1]
        density = inputs[2]
        allLines = inputs[3]
        thisIon = ch.ion(ionS, temperature, density)
        thisIon.intensity(allLines = allLines)
        outList = [ionS, thisIon.Intensity]
        outQueue.put(outList)
    return
    #
    #-----------------------------------------------------
    #
def makeMatchPkl(specData, temp, dens, wghtFactor = 0.25,  abundanceName = None, minAbund=1.e-6, useMgofnt=1, verbose=0):
    '''
    input a data dictionary and instantiate a dem class,
    and run mgofnt and then make a pickle file
    to use multiprocessing, this needs to be run in an ipython console

    Parameters

    specData : dict
        the observed line intensities, wavelegths ...
    '''

    mydem = maker(specData, wghtFactor = wghtFactor, abundanceName =  abundanceName, minAbund=minAbund, verbose=verbose)
    if useMgofnt:
        mydem.mgofnt(temp, dens, verbose=1)
    else:
        mydem.gofnt(temp, dens)
#    mydem.predictor()
    filename = specData['filename']
    outname = os.path.splitext(filename)[0] + '-' + mydem.AbundanceName +'.pkl'
    mydem.PklName = outname
    print(' pickle name = %s'%(outname))
    outpt = open(outname, 'wb')
    pickle.dump(mydem, outpt)
    outpt.close()
    print(' should now move this to the working dir where it will be used')
    #
    #-----------------------------------------------------
    #
class maker(ionTrails):
    '''
    a class matching observed lines to lines in the CHIANTI database
    '''
    def __init__(self, specData, wghtFactor = 0., ionList=False, allLines=True, abundanceName = None, minAbund=1.e-6, verbose=False):
        '''
        input a list of wavelengths and a wavelength difference
        find a list of predicted spectral lines for each wavelength
        wgtFactor:  multiply this by the predicted to get the weight of the intensity difference
            applied in getWeighted()
        all = 0 -> only previously observed wavelengths are used
        dwvl = wavelength difference between CHIANTI and the observed lines
            if == 0 then use use that in the observed data
        Parameters
        ----------
        specData : dict
            contains the following keys
            intensity - a list of observed line intensities
            wvlObs - a list of observed wavelengths
            identities, dwvl
        wghtFactor : float
            the factor used in calculating the weight chi-squared
        ionList :  list
            a list of ions, such as fe_14, to be searched
        allLines : bool
            if true, unobserved lines in CHIANTI are included in the search
        abundanceName : str
            the name of the elemental abundance file to be used,
            if not set, the default abundance file is used
        minAbund : float
            sets the minimum abundance for an element to be included in the search
        verbose : bool
            if True, additional output is sent to the terminal
        '''
        #
        # --------------------------------------------------------------------------------
        #
        self.Defaults = chdata.Defaults
        self.XUVTOP = os.environ['XUVTOP']
        print(' XUVTOP = %s'%(self.XUVTOP))
        if abundanceName is not None:
            if abundanceName in chdata.AbundanceList:
                self.AbundanceName = abundanceName
            else:
                abundChoices = chdata.AbundanceList
                abundChoice = chGui.gui.selectorDialog(abundChoices,label='Select Abundance name', multiChoice=False)
                abundChoice_idx = abundChoice.selectedIndex
                self.AbundanceName = abundChoices[abundChoice_idx[0]]
        else:
            self.AbundanceName = self.Defaults['abundfile']
        print(' abundanceName = %s'%(self.AbundanceName))

        self.SpecData = specData
        self.MinAbund = minAbund
        print(' minimum abundance = %10.2e'%(minAbund))
        mlInfo = io.masterListInfo()
        self.Intensity = specData['intensity']
        self.IonS = specData['ions']
        self.IonSet = set(specData['ions'])
        self.Dwvl = specData['dwvl']
        reduceNobs = 0.
        for anion in self.IonSet:
            same = specData['ions'].count(anion)
            reduceNobs += np.sqrt(float(same))
        self.ReduceNobs = reduceNobs
        self.Nobs = len(specData['wvl0'])
        print(' # of observables / reduce # = %i  %10.2f'%(self.Nobs,self.ReduceNobs))
        self.Wvl = specData['wvl0']
        #
        self.SpecData['wghtFactor'] = []
        if wghtFactor > 0.:
            for iwvl in range(self.Nobs):
                self.WghtFactor = self.SpecData['intStd'][iwvl]/self.SpecData['intensity'][iwvl] + wghtFactor
        else:
            for iwvl in range(self.Nobs):
                self.WghtFactor = self.SpecData['intStd'][iwvl]/self.SpecData['intensity'][iwvl]
        #
        self.AllLines = allLines
        if 'intensity' in specData.keys():
            self.Intensity = specData['intensity']
        masterlist = io.masterListRead()
#        matches = [{'ion':[], 'wvl':[]}]*len(wvl) - it won't work this way
        matches = []
        for iwvl in range(len(self.Wvl)):
            matches.append({'ion':[], 'wvl':[], 'lineIdx':[], 'wvldiff':[], 'lvl1':[], 'lvl2':[], 'pretty1':[], 'pretty2':[], 'predictedLine':[],'iPredictedLine':0, 'obsIntensity':self.Intensity[iwvl], 'exptIon':self.IonS[iwvl], 'obsWvl':self.Wvl[iwvl]})
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
                        print('')
            masterlist = alist
        #
        if minAbund:
            self.MinAbund = minAbund
            abundanceName = chdata.Defaults['abundfile']
            abundanceAll = io.abundanceRead(abundancename = abundanceName)
            revList = []
            for one in masterlist:
                params = util.convertName(one)
                z = params['Z']
                if abundanceAll['abundance'][z-1] > minAbund:
                    revList.append(one)
            masterlist = revList
        for thision in masterlist:
            if verbose:
                print(' - - - - - - - - - - - ')
                print(' thision = %s'%(thision))
            btw = util.between(self.Wvl, [mlInfo[thision]['wmin'], mlInfo[thision]['wmax']])
            if len(btw):
                wgfa = io.wgfaRead(thision)
                cnt = wgfa['wvl'].count(0.)
                # in order to match with emiss
                for icnt in range(cnt):
                    idx = wgfa['wvl'].index(0.)
                    wgfa['wvl'].pop(idx)
                    wgfa['lvl1'].pop(idx)
                    wgfa['lvl2'].pop(idx)
                    wgfa['pretty1'].pop(idx)
                    wgfa['pretty2'].pop(idx)
                thesewvl = np.asarray(wgfa['wvl'])
                nonzed = thesewvl != 0.
                thesewvl =thesewvl[nonzed]
                if allLines:
                    thesewvl = np.abs(thesewvl)
    #            matches[iwvl]['lvl1'] = []
    #            matches[iwvl]['lvl2'] = []
    #            matches[iwvl]['pretty1'] = []
    #            matches[iwvl]['pretty2'] = []
                for iwvl, awvl in enumerate(self.Wvl):
#                    if verbose:
#                        print(' iwvl,awvl = %5i  %10.3f'%(iwvl, awvl))
#                        print('  in matches[iwvl][ion] = %s'%(matches[iwvl]['ion']))
#                        print('  in matches[iwvl][wvl] = %10.3f'%(matches[iwvl]['wvl']))
                    diff = np.abs(thesewvl - awvl)
                    good = diff < self.Dwvl[iwvl]
                    ngood = good.sum()
                    if ngood:
    #                    print '    ion = ', thision
    #                    print '    wvl = ', thesewvl[good].tolist()
                        matches[iwvl]['ion'].append(thision)
                        matches[iwvl]['wvl'].append(thesewvl[good].tolist())
                        matches[iwvl]['lineIdx'].append(np.arange(len(thesewvl))[good].tolist())
                        matches[iwvl]['wvldiff'].append(diff[good])
                        alvl1 = []
                        alvl2 = []
                        ap1 = []
                        ap2 = []
                        for adx in matches[iwvl]['lineIdx'][-1]:
    #                        print(' iwvl = %5i adx = %5i'%(iwvl, adx))
                            alvl1.append(wgfa['lvl1'][adx])
                            alvl2.append(wgfa['lvl2'][adx])
                            ap1.append(wgfa['pretty1'][adx])
                            ap2.append(wgfa['pretty2'][adx])
                        matches[iwvl]['lvl1'].append(alvl1)
                        matches[iwvl]['lvl2'].append(alvl2)
                        matches[iwvl]['pretty1'].append(ap1)
                        matches[iwvl]['pretty2'].append(ap2)
#                        if verbose:
#                            print('     out matches[iwvl][ion] = %s'%(matches[iwvl]['ion']))
    #                    print '     out matches[iwvl][wvl] = ', matches[iwvl]['wvl']
    #                    print ' good ion lines = ', thesewvl[good]
    #                print ' thisMatch = ',  thisMatch
    #                matches[iwvl].append(thisMatch)
            else:
                if verbose:
                    print(' no lines in this wavelength range')
        for amatch in matches:
            nions = len(amatch['ion'])
            amatch['predictedLine'] = [0]*nions
        self.match = matches

        #em = None,
    def argCheck(self, temperature=None, eDensity=None, pDensity='default', verbose=0):
        ''' to check the compatibility of the three arguments
        and put them into numpy arrays of atleast_1d
        '''
        if temperature is not None:
            self.Temperature = np.atleast_1d(temperature)
            if isinstance(self.Temperature[0], str):
                raise ValueError(' temperature can not be a string')
            if np.any(self.Temperature <= 0.):
                raise ValueError(' all temperatures must be positive')
            self.Ntemp = self.Temperature.size
        else:
            raise ValueError('temperature not defined')

        if pDensity is 'default':
            self.p2eRatio()

        if eDensity is not None:
            self.EDensity = np.atleast_1d(eDensity)
            if isinstance(self.EDensity[0], str):
                raise ValueError(' EDensity can not be a string')
            if np.any(self.EDensity <= 0.):
                raise ValueError(' all densities must be positive')
            self.Ndens = self.EDensity.size
        # needed when doing ioneq.calculate()
        else:
            self.Ndens = 0

        if verbose:
            print('Temperature size:  %5i  Density size:  %5i'%(self.Ntemp, self.Ndens))

        self.NTempDens = max(self.Ndens,self.Ntemp)
        if self.Ndens > 1 and self.Ntemp == 1:
            self.Temperature = np.tile(self.Temperature, self.NTempDens)
        elif self.Ndens == 1 and self.Ntemp > 1:
            self.EDensity = np.tile(self.EDensity, self.NTempDens)

        if hasattr(self,'EDensity') and hasattr(self,'Temperature') and self.Temperature.size != self.EDensity.size:
            raise ValueError('Temperature and density must be the same size.')
        if pDensity is not None:
            if pDensity is 'default' and eDensity is not None:
                self.PDensity = self.ProtonDensityRatio*self.EDensity
            else:
                self.PDensity = np.atleast_1d(pDensity)
                if self.PDensity.size < self.Ndens:
                    np.tile(self.PDensity, self.Ndens)
                    self.NpDens = self.NpDens.size

        #
        # -----------------------------------------------------------------
        #
    def gofnt(self, temperature, density, verbose=1):
        '''
        calculate the gofnt function for each of the matched lines
        do each ion only once
        '''
        t1 = datetime.now()
        self.XUVTOP = os.environ['XUVTOP']
        #  quick and dirty
        #   em = None,
        self.argCheck(temperature=temperature, eDensity=density, pDensity=None,  verbose=0)
        temperature = self.Temperature
        density = self.EDensity

#        if not np.iterable(temperature):
#            temperature = [temperature]
#        if not np.iterable(density):
#            density = [density]
#        self.Temperature = np.asarray(temperature)
#        self.Density = np.asarray(density)
#        nTempDens = max([len(temperature), len(density)])
#        if nTempDen == 1:
#            print(' the number of temperatures or densities should be greater than 1')
#            return
#        self.NTempDen = nTempDen
        nTempDens = self.NTempDens
        if verbose:
            print(' temperature size:  %5i'%(self.Temperature.size))
            print(' density     size:  %5i'%(self.EDensity.size))
        ionList = []
        for iwvl, amatch in enumerate(self.match):
            for someIon in amatch['ion']:
#                print 'someIon = ', someIon
                if someIon not in ionList:
                    ionList.append(someIon)
        for iwvl, amatch in enumerate(self.match):
            nPredictedLines = 0
            for awvl in amatch['wvl']:
                nPredictedLines +=  len(awvl)
            self.match[iwvl]['intensity'] = np.zeros((nPredictedLines,  nTempDens), 'float64')

        self.ionList = ionList
        for iwvl in range(len(self.match)):
            self.match[iwvl]['intensitySum'] = np.zeros(nTempDens, 'float64')
        for someIon in ionList:
            # already know this ion is needed
            #if verbose:
                #print(' someIon = %s'%(someIon))
            thisIon = ch.ion(someIon, temperature, density, abundance = self.AbundanceName)
            if hasattr(thisIon, 'IoneqOne'):
                thisIon.intensity()
                intensity = thisIon.Intensity['intensity']
                #
#                if len(intensity.shape) == 1:
#                    nTempDens = 1
#                else:
#                    nTempDens = intensity.shape[0]
    #            print(' nTempDens = ', nTempDens)
                for iwvl, amatch in enumerate(self.match):
    #                self.match[iwvl]['intensitySum'] = np.zeros(nTempDens, 'float64')
                    #  this is data for each line
                    if someIon in amatch['ion']:
                        kon = amatch['ion'].index(someIon)
                        if verbose:
                            print(' using %s'%(someIon))
    #                for kon, anIon in enumerate(amatch['ion']):
    #
        #                print ' linIdx = ', amatch['lineIdx'], amatch
                        predictedLine = []
                        for aline in amatch['lineIdx'][kon]:
                            #print(' %s   %6i  %12.3f '%( someIon, aline, thisIon.Intensity['wvl'][aline]))
                            self.match[iwvl]['intensitySum'] += intensity[:, aline]
                            iPredictedLine = self.match[iwvl]['iPredictedLine']
                            predictedLine.append(iPredictedLine)
                            self.match[iwvl]['intensity'][iPredictedLine] = intensity[:, aline]
                            self.match[iwvl]['iPredictedLine'] += 1
                        self.match[iwvl]['predictedLine'][kon] = predictedLine
                self.Tmax = np.zeros_like(self.Wvl)
#                self.Dmax = np.zeros_like(self.Wvl)
                for iwvl, amatch in enumerate(self.match):
                    gfun = amatch['intensitySum']
                    peak = gfun == gfun.max()
#                    if len(temperature) > 1:
                    self.match[iwvl]['tmax'] = temperature[peak]
#                    elif len(density) >1:
#                        self.match[iwvl]['dmax'] = density[peak]
            else:
               print(' IoneqOne not available for %s'%(someIon))
        t2 = datetime.now()
        dt = t2 - t1
        print(' elapsed seconds = %12.3f'%(dt.seconds))
        #for iwvl, amatch in enumerate(self.match):
            #gfun = amatch['intensitySum']
            #peak = gfun == gfun.max()
            #self.match[iwvl]['tmax'] = temperature[peak]
        #
        # ---------------------------------------------------------------------
        #
    def mgofnt(self, temperature, density, proc=6,  timeout=0.1, verbose=0):
        '''
        calculate the gofnt function for each of the matched lines
        this is the multiprocessing version
        do each ion only once
        '''
        t1 = datetime.now()
        self.XUVTOP = os.environ['XUVTOP']
        #  quick and dirty
        self.argCheck(temperature=temperature, eDensity=density, pDensity=None,  verbose=0)

        temperature = self.Temperature
        density = self.EDensity
        nTempDens = self.NTempDens

        for iwvl, amatch in enumerate(self.match):
            nPredictedLines = 0
            for awvl in amatch['wvl']:
                nPredictedLines +=  len(awvl)
            self.match[iwvl]['intensity'] = np.zeros((nPredictedLines,  nTempDens), 'float64')

        ionList = []
        for iwvl, amatch in enumerate(self.match):
            for someIon in amatch['ion']:
#                print 'someIon = ', someIon
                if someIon not in ionList:
                    ionList.append(someIon)

        self.ionList =ionList

        for iwvl in range(len(self.match)):
            self.match[iwvl]['intensitySum'] = np.zeros(nTempDens, 'float64')
        #
        #  ion multiprocessing setup
        #
        ionWorkerQ = mp.Queue()
        ionDoneQ = mp.Queue()
        proc = min([proc, mp.cpu_count()])
        #
        #
        self.Todo = []
        for someIon in ionList:
            # already know this ion is needed
            print(' someIon = %s'%(someIon))
            ionWorkerQ.put((someIon, temperature, density, self.AllLines))
            self.Todo.append(someIon)
        #
        ionWorkerQSize = ionWorkerQ.qsize()
        ionProcesses = []
        if ionWorkerQSize < proc:
            nproc = ionWorkerQSize
        for i in range(proc):
            p = mp.Process(target=doDemGofntQ, args=(ionWorkerQ, ionDoneQ))
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
        self.Finished = []
        #
        for ijk in range(ionWorkerQSize):
            out = ionDoneQ.get()
            someIon = out[0]
            print('processing ion = %s'%(someIon))
            self.Finished.append(someIon)
            intensity = out[1]['intensity']

            for iwvl, amatch in enumerate(self.match):
#                self.match[iwvl]['intensitySum'] = np.zeros(nTempDens, 'float64')
                #  this is data for each line
                if someIon in amatch['ion']:
                    kon = amatch['ion'].index(someIon)
#                    if verbose:
#                        print('kon: %5i using %s'%(kon, someIon))
#                        if 'errorMessage' in out[1].keys():
#                            print(' in mgofnt %s'%(out[1]['errorMessage']))

                    predictedLine = []
                    for aline in amatch['lineIdx'][kon]:
    #                    print ' ion, lineIdx = ', anIon, aline
                        self.match[iwvl]['intensitySum'] += intensity[:, aline]
                        iPredictedLine = self.match[iwvl]['iPredictedLine']
                        predictedLine.append(iPredictedLine)
                        self.match[iwvl]['intensity'][iPredictedLine] = intensity[:, aline]
                        self.match[iwvl]['iPredictedLine'] += 1
                    self.match[iwvl]['predictedLine'][kon]=predictedLine
        #
        self.Tmax = np.zeros_like(self.Wvl)
        nT = self.Temperature.size
        trange = np.arange(nT)
        minIdx = []
        maxIdx = []
        for iwvl, amatch in enumerate(self.match):
            gfun = amatch['intensitySum']
            peak = gfun == gfun.max()
            if len(temperature) > 1:
                self.match[iwvl]['tmax'] = temperature[peak]
            elif len(density) >1:
                self.match[iwvl]['dmax'] = density[peak]
#            self.Tmax[iwvl] = temperature[peak]
            good = amatch['intensitySum'] > 0.
            if good.sum() > 0.:
                tRangeGood = trange[good]
                minIdx.append(tRangeGood.argmin())
                maxIdx.append(tRangeGood.argmax())
            else:
                print(' intensitySum == 0. for iwvl %5i'%(iwvl))
        self.MinIndex = max(minIdx)
        self.MaxIndex = min(maxIdx)


        #
        #
        for p in ionProcesses:
            if not isinstance(p, str):
                p.terminate()

        t2 = datetime.now()
        dt = t2 - t1
        print(' elapsed seconds = %12.3f'%(dt.seconds))
        #
        # --------------------------------------------------------
        #
    def emNtSetIndices(self, indices, verbose=0):
        '''
        to set the indices of the N temperature EM distribution
        '''
        print(' should use emSetIndices')
        if hasattr(self, 'Temperature'):
            self.EmIndices = indices
            if verbose:
                for adx in indices:
                    print('index %5i temperature %12.2e'%(adx, self.Temperature[adx]))
        else:
            print(' need to call gofnt/mgofnt first')
        #
        # --------------------------------------------------------
        #
    def emSetIndices(self, indices, verbose=0):
        '''
        to set the indices of the N temperature/density EM distribution
        '''
        if hasattr(self, 'Temperature'):
            self.EmIndices = np.atleast_1d(indices)
            if verbose:
                for adx in self.EmIndices:
                    print('index %5i temperature %12.2e'%(adx, self.Temperature[adx]))
        else:
            print(' need to call gofnt/mgofnt first')
        #
        #------------------------------------------
        #
    def emSet(self, value):
        '''
        sets the EM values for a N temperature EM distribution
        '''
        emValue = np.atleast_1d(value)
        if hasattr(self, 'EmIndices') and len(emValue) == len(self.EmIndices):
            em = np.zeros_like(self.Temperature)
            for i, adx in enumerate(self.EmIndices):
                em[adx] = 10.**emValue[i]
            self.Em = em
        else:
            print('in emNtexp, either EmIndices not set or no. of values != no. of indices ')
            if hasattr(self, 'EmIndices'):
                print(' # of indices = %i, len(EmIndices = %i'%(len(value),len(self.EmIndices)))
            else:
                print(' self.EmIndices does not exist')
                print('# of indices = '%(len(value)))
        #
        # ---------------------------------------------------------
        #
    def emFitPlot(self):
        '''
        to plot the emission measures derived from search over temperature
        '''
        if not hasattr(self, 'SearchData'):
            print(' must run search*t... first')
            return
        else:
            em = self.SearchData['best']['em']
            temp = self.SearchData['best']['temperature']
#            for it, em1 in enumerate(em):
#                plt.loglog([temp[it], temp[it]], [em1/1000., em1], '-k', linewidth=2)
            plt.loglog([temp, temp], [em/1000., em], '-k', linewidth=2)
            plt.loglog([float(self.Temperature.min()), float(self.Temperature.max())], [em, em], '-k', linewidth=2 )
        #
        # ---------------------------------------------------------
        #
    def emPlot(self, vs='T'):
        '''
        to plot line intensities divided by gofnt
        '''
        nInt = len(self.Intensity)
#        print(' nInt = %5i'%(nInt))
        if not hasattr(self, 'Temperature'):
            print(' must run mgofnt or gofnt first')
            return
        elif vs is 'T':
            if self.Ntemp > 1:
                em = np.zeros((nInt, self.Ntemp), 'float64')
                xvar = self.Temperature
            for idx in range(nInt):
                nonzed = self.match[idx]['intensitySum'] > 0.
                large = self.match[idx]['intensitySum'] >  self.match[idx]['intensitySum'].max()*0.01
                realgood = np.logical_and(nonzed, large)
    #            print(' for idx= %5i size of realgood = %5i'%(idx, realgood.size))
                if realgood.sum() > 0 and self.Intensity[idx] > 0.:
                    em[idx][realgood] = self.Intensity[idx]/self.match[idx]['intensitySum'][realgood]
                    plt.loglog(xvar[realgood], em[idx][realgood], lw=2)
        #            minv = em[idx][nonzed] == em[idx][nonzed].min()
                    if len(self.Temperature) > 1:
                        lblx = self.match[idx]['tmax']
                    else:
                        lblx = self.match[idx]['dmax']
        #            lbly = em[idx][nonzed][minv[0]]
                    lbly = np.min(em[idx][realgood])*1.05
                    wvlstr = ' %8.3f'%(self.Wvl[idx])
        ##            lbl = self.match[idx]['ion'] + wvlstr
    #                print ' lblx,lbly,wvlstr, minv = ', lblx, lbly, wvlstr
                    plt.text(lblx, lbly, wvlstr.strip())
                else:
                    print('no values for idx = %5i wvl = %8.2f'%(idx, self.Wvl[idx]))
                    print(' nonzed = %i'%(nonzed.sum()))
                    print('intensity = %10.2e'%(self.Intensity[idx]))
            plt.ylabel('Emission Measure (cm$^{-5}$)', fontsize=14)
            plt.xlabel('Temperature (K)',  fontsize=14)
            if 'file' in self.SpecData.keys():
                plt.title(self.SpecData['file'])
        elif vs is not 'T':
            if self.Ndens > 1:
                em = np.zeros((nInt, self.Ndens), 'float64')
                xvar = self.EDensity
            for idx in range(nInt):
                nonzed = self.match[idx]['intensitySum'] > 0.
                large = self.match[idx]['intensitySum'] >  self.match[idx]['intensitySum'].max()*0.01
                realgood = np.logical_and(nonzed, large)
    #            print(' for idx= %5i size of realgood = %5i'%(idx, realgood.size))
                if realgood.sum() > 0 and self.Intensity[idx] > 0.:
                    em[idx][realgood] = self.Intensity[idx]/self.match[idx]['intensitySum'][realgood]
                    plt.loglog(xvar[realgood], em[idx][realgood], lw=2)
        #            minv = em[idx][nonzed] == em[idx][nonzed].min()
#                    lblx = self.match[idx]['dmax']
                    lblx = self.EDensity
        #            lbly = em[idx][nonzed][minv[0]]
#                    lbly = np.min(em[idx][realgood])
                    lbly = em[idx][realgood]*1.05
                    wvlstr = ' %8.3f'%(self.Wvl[idx])
        ##            lbl = self.match[idx]['ion'] + wvlstr
    #                print ' lblx,lbly,wvlstr, minv = ', lblx, lbly, wvlstr
                    plt.text(lblx[0], lbly[0], wvlstr.strip())
                    plt.text(lblx[-1], lbly[-1], wvlstr.strip(), horizontalalignment='right')
                else:
                    print('no values for idx = %5i wvl = %8.2f'%(idx, self.Wvl[idx]))
                    print(' nonzed = %i'%(nonzed.sum()))
                    print('intensity = %10.2e'%(self.Intensity[idx]))
            plt.ylabel('Emission Measure (cm$^{-5}$)', fontsize=14)
            plt.xlabel('Electron Density (cm$^{-3}$)',  fontsize=14)
            if 'file' in self.SpecData.keys():
                plt.title(self.SpecData['file'])
        #
        # --------------------------------------------------------------------------
        #
    def diffPrintChi(self, dir = '.', filename='diffPrintChi.txt'):
        '''
        calculates the weighted and straight differences between observed and predicted
        prints the values saves the as a dictionary self.Diff
        to be used together with a prior brute-force chi-squared minimization approach
        '''
        wDiff = []
        relDevList = []
        stdDevList = []
        intOverPred = []
        dash = ' -------------------------------------------------'
        pformat1 = ' %5i %7s %10.3f %10.2e %10.2e %10.3f %10.3f %10.3f %10.3f'
        pformat1o = ' %5i %7s %10.3f %10.2e %10.2e %10.3f %10.3f %10.3f \n'
        sformat  = ' %5s %7s %10s %10s %10s %10s %10s %10s %10s'
        with open(filename, 'w') as outpt:

            try:
                idx = self.SearchData['best']['idx']
                dens = self.SearchData['best']['density']
                temp = self.SearchData['best']['temperature']
                emfit = self.SearchData['best']['emfit']
                em = self.SearchData['best']['em']
                print(' %5i %10.2e %10.2e %10.3f %10.2e'%(idx, dens, temp, emfit, em))
                outpt.write(' %5i %10.2e %10.2e %10.3f %10.2e \n'%(idx, dens, temp, emfit, em))
            except:
                emIndices = self.EmIndices
                for emidx in emIndices:
                    print(' %5i %12.2e %12.2e %12.2e'%(emidx, self.EDensity[emidx],  self.Temperature[emidx], self.Em[emidx]))
                    outpt.write(' %5i %12.2e %12.2e %12.2e \n'%(emidx, self.EDensity[emidx],  self.Temperature[emidx], self.Em[emidx]))
            finally:
                pass
            print(dash)
            outpt.write(dash +'\n')
            print(' chi = abs(int/(2*wght*pred))  relDiff = (int - pred)/pred')
            outpt.write(' chi = abs(int/(2*wght*pred))  relDiff = (int - pred)/pred \n')
            sprint = sformat%('iwvl', 'ionS', 'wvl', 'intensity', 'predicted', 'int/pred', 'chi', 'relDiff', 'stdDev')
            print(sprint)
            outpt.write(sprint +'\n')
            for iwvl,  amatch in enumerate(self.match):
                if amatch['predicted'] > 0. :
                    chi = np.abs(self.Intensity[iwvl]/(2.*self.WghtFactor*amatch['predicted']))
                    wDiff.append(chi)
                    intOverPred.append(np.abs(self.Intensity[iwvl]/amatch['predicted'] - 1.))
                    relDevList.append(np.abs(self.Intensity[iwvl]-amatch['predicted'])/amatch['predicted'])
                    stdDevList.append(relDevList[-1]**2)
                    pstring = pformat1%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'], self.Intensity[iwvl]/amatch['predicted'], chi, relDevList[-1], stdDevList[-1] )
                else:
                    pstring = pformat1%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'], -1.)
                print(pstring)
                outpt.write(pstring + '\n')
            relDev = np.asarray(relDevList, np.float64)
            stdDev = np.asarray(stdDevList, np.float64)
            print(' mean of relative difference = %10.3f std = %10.3f stdDev  %10.3f'%(relDev.mean(), relDev.std(), stdDev.mean()))
            outpt.write(' mean of relative difference = %10.3f std = %10.3f stdDev  %10.3f \n'%(relDev.mean(), relDev.std(), stdDev.mean()))
            print(' mean of intOverPred = abs(int/pred -1.) %10.3f'%(np.mean(intOverPred)))
            outpt.write(' mean of intOverPred = abs(int/pred -1.) %10.3f \n'%(np.mean(intOverPred)))
            threeSig = 3.*relDev.std()
            print(' 3*std = %10.3f'%(threeSig))
            poor = relDev > threeSig
            idx = np.arange(relDev.size)
            pdx = idx[poor]
            for i in pdx:
                print('%5i %s %10.3f %10.3f %10.3f'%(i, self.IonS[i],  self.Wvl[i], relDev[i], stdDev[i]))
                outpt.write('%5i %s %10.3f %10.3f %10.3f \n'%(i, self.IonS[i],  self.Wvl[i], relDev[i], stdDev[i]))
        diffDict = {'diff':relDev, 'wvl':self.Wvl, 'ionS':self.IonS, '3sig':threeSig, 'poor':poor, 'stdDev':stdDev}
        self.SearchData['diff'] = diffDict
        #
        # --------------------------------------------------------------------------
        #
    def diffPrintMc(self, dir = '.', filename='diffPrintMc.txt'):
        '''
        calculates the weighted and straight differences between observed and predicted
        prints the values saves to a file from a PyMC run
        '''
        wDiff = []
        relDevList = []
        stdDevList = []
        intOverPred = []
        dash = ' -------------------------------------------------'
        pformat1 = ' %5i %7s %10.3f %10.2e %10.2e %10.3f %10.3f %10.3f %10.3f'
        pformat1a = ' %5i %7s %10.3f %10.2e %10.2e *****NID******'
        sformat  = ' %5s %7s %10s %10s %10s %10s %10s %10s %10s'
        fsformat  = ' %5s %7s %10s %10s %10s %10s %10s %10s %10s\n'
        with open(filename, 'w') as outpt:
            emIndices = self.EmIndices
            print(' %5s %12s %12s %12s'%('index',  'density',  'temperature',  'Em'''))
            outpt.write(' %5s %12s %12s %12s \n'%('index',  'density',  'temperature',  'Em'''))
            for emidx in emIndices:
                print(' %5i %12.2e %12.2e %12.2e'%(emidx, self.EDensity[emidx],  self.Temperature[emidx], self.Em[emidx]))
                outpt.write(' %5i %12.2e %12.2e %12.2e \n'%(emidx, self.EDensity[emidx],  self.Temperature[emidx], self.Em[emidx]))
            print(dash)
            outpt.write(dash+'\n')
            print(' chi = abs(int/(2*wght*pred))  strDiff = (int - pred)/pred')
            outpt.write(' chi = abs(int/(2*wght*pred))  strDiff = (int - pred)/pred \n')
            print(sformat%('iwvl', 'ionS', 'wvl', 'intensity', 'predicted', 'int/pred', 'chi', 'relDev', 'stdDev'))
            outpt.write(fsformat%('iwvl', 'ionS', 'wvl', 'intensity', 'predicted', 'int/pred', 'chi', 'relDev',  'stdDev'))
            for iwvl,  amatch in enumerate(self.match):
                if amatch['predicted'] > 0. :
                    chi = np.abs(self.Intensity[iwvl]/(2.*self.WghtFactor*amatch['predicted']))
                    wDiff.append(chi)
                    intOverPred.append(np.abs(self.Intensity[iwvl]/amatch['predicted'] - 1.))
                    relDevList.append(np.abs(self.Intensity[iwvl]-amatch['predicted'])/self.Intensity[iwvl])
                    stdDevList.append(relDevList[-1])
                    pstring = pformat1%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'], self.Intensity[iwvl]/amatch['predicted'], chi, relDevList[-1],  stdDevList[-1] )
                else:
                    pstring = pformat1a%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'])
                print(pstring)
                outpt.write(pstring+'\n')
            diff = np.asarray(relDevList, np.float64)
            stdDev = np.asarray(stdDevList, np.float64)
            print(' mean of relative Deviation = %10.3f 3*relDev = %10.3f stdDev = %10.3f\n'%(diff.mean(), 3.*diff.mean(), stdDev.mean()))
            outpt.write(' mean of rel Dev = %10.3f 3*rel Dev  = %10.3f stdDiff:  %10.3f \n'%(diff.mean(), 3.*diff.mean(), stdDev.mean()))
            print(' mean of intOverPred = abs(int/pred -1.) %10.3f'%(np.mean(intOverPred)))
            outpt.write(' mean of intOverPred = abs(int/pred -1.) %10.3f \n'%(np.mean(intOverPred)))
            threeSig = 3.*diff.mean()
            print(' 3*std = %10.3f'%(threeSig))
            outpt.write(' 3*std = %10.3f \n'%(threeSig))
            poor = diff > threeSig
            idx = np.arange(diff.size)
            pdx = idx[poor]
            for i in pdx:
                print('%5i %s %10.3f %10.3f'%(i, self.IonS[i],  self.Wvl[i], diff[i]))
                outpt.write('%5i %s %10.3f %10.3f \n'%(i, self.IonS[i],  self.Wvl[i], diff[i]))
        self.diffDict = {'diff':diff, 'wvl':self.Wvl, 'ionS':self.IonS, '3sig':threeSig, 'poor':poor}
#        self.SearchData['diff'] = diffDict
        #
        # --------------------------------------------------------------------------
        #
    def predict(self):
        '''
        to predict the intensities of the observed lines from an emission measure
        the emission measure is already specified as self.Em which is an np array
        '''
        #
        for iwvl, amatch in enumerate(self.match):
            pred = np.sum(amatch['intensitySum']*self.Em)
            try:
                self.match[iwvl]['predicted'] = pred
            except:
                self.match[iwvl]['predicted'] = 0.
                print(' error in predict, iwvl = %5i'%(iwvl))
        #
        # --------------------------------------------------------------------------
        #
    def predictPrint(self, minContribution=0.1, outfile=0, verbose=0):
        '''
        to predict the intensities of the observed lines from an emission measure
        the emission measure is already specified as self.Em which is an np array
        '''
        dash = ' -------------------------------------------------'
        pformat1 = ' %5i %7s %10.3f %10.2e %10.2e %10.3f %10.3f'
        pformat1s = ' %5s %7s %10s %10s %10s %10s %10s'
        pformat2 = '         %s'
        pformat3 = '        %10.3f %4i %4i %20s - %20s %5i %5i %7.3f'
        pformat3s = '        %10s %4s %4s %20s - %20s %5s %5s %s'
        try:
            emIndices = self.EmIndices
            for emidx in emIndices:
                print(' %5i %12.2e %12.2e %12.2e'%(emidx, self.EDensity[emidx],  self.Temperature[emidx], self.Em[emidx]))
        except:
            print(' attribute EmIndices has not been created')
            return
        print(dash)
        pstring1 = pformat1s%('iwvl', 'IonS', 'wvl', 'Int',  'Pred', 'Int/Pred', 'chi')
        print(pstring1)
        pstring3 = pformat3s%('wvl', 'lvl1', 'lvl2', 'lower', 'upper', 'lineIdx',  'predLine', 'contribution')
        print(pstring3)
        print(dash)
        for iwvl,  amatch in enumerate(self.match):
            if amatch['predicted'] > 0. :
                chi = np.abs(self.Intensity[iwvl]/(2.*self.WghtFactor*amatch['predicted']))
                pstring = pformat1%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'], self.Intensity[iwvl]/amatch['predicted'], chi  )
            else:

                pstring = pformat1%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'], -1.)
            print(pstring)
            #
            # now check line contributions
            #
            for jon,  anion in enumerate(amatch['ion']):
                ionPrint = 0
                for iline, awvl in enumerate(amatch['wvl'][jon]):
                    contrib = (amatch['intensity'][amatch['predictedLine'][jon][iline]]*self.Em).sum()/amatch['predicted']
                    if contrib > minContribution:
                        if ionPrint == 0:
                            print(pformat2%(anion))
                            ionPrint = 1
                        print(pformat3%(awvl, amatch['lvl1'][jon][iline], amatch['lvl2'][jon][iline], amatch['pretty1'][jon][iline], amatch['pretty2'][jon][iline].ljust(20), amatch['lineIdx'][jon][iline], amatch['predictedLine'][jon][iline], contrib))
                        print(dash)
            print(dash)
        if outfile:
            with open(outfile, 'w') as outpt:
                try:
                    idx = self.SearchData['best']['idx']
                    dens = self.SearchData['best']['density']
                    temp = self.SearchData['best']['temperature']
                    emfit = self.SearchData['best']['emfit']
                    em = self.SearchData['best']['em']
                    print(' %5i %10.2e %10.2e %10.3f %10.2e'%(idx, dens, temp, emfit, em))
                except:
                    for emidx in self.EmIndices:
                        outpt.write(' %5i%12.2e %12.2e \n'%(emidx, self.Temperature[emidx], self.Em[emidx]))
                finally:
                    pass
                outpt.write(dash +'\n')
                pstring1 = pformat1s%('iwvl', 'IonS', 'wvl', 'Int',  'Pred', 'Int/Pred', 'chi')
                outpt.write(pstring1 +'\n')
                outpt.write(pstring3 + '\n')
                for iwvl,  amatch in enumerate(self.match):
                    chi = np.abs(self.Intensity[iwvl]/(2.*self.WghtFactor*amatch['predicted']))
                    pstring = pformat1%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'], self.Intensity[iwvl]/amatch['predicted'], chi  )
                    outpt.write(pstring +'\n')
                    #
                    # now check line contributions
                    #
                    for jon,  anion in enumerate(amatch['ion']):
                        ionPrint = 0
                        for iline, awvl in enumerate(amatch['wvl'][jon]):
                            contrib = (amatch['intensity'][amatch['predictedLine'][jon][iline]]*self.Em).sum()/amatch['predicted']
                            if contrib > minContribution:
                                if ionPrint == 0:
                                    outpt.write(pformat2%(anion)+'\n')
                                    ionPrint = 1
                                outpt.write(pformat3%(awvl, amatch['lvl1'][jon][iline], amatch['lvl2'][jon][iline], amatch['pretty1'][jon][iline], amatch['pretty2'][jon][iline], amatch['lineIdx'][jon][iline], amatch['predictedLine'][jon][iline], contrib) +'\n')
                                outpt.write(dash +'\n')
                    outpt.write(dash +'\n')

            # now with latex notation

                outpt.write(dash +'\n')
                outpt.write(' ------------------- LATEX ---------------------\n')
                outpt.write(dash +'\n')
                pformat1S = ' %5s %10s & %10s & %10s & %10s & %10s  \n'
                pstring1 = pformat1S%('iwvl', 'wvl', 'wvl', 'Int',  'Pred', 'Int/Pred' )
                outpt.write(pstring1 +'\n')
                outpt.write(dash +'\n')
                pformat1L = ' %5i %10.3f & %10.3f & %10.2e & %10.2e & %10.3f  \n'

                pformat3L = '        %10.3f & %10.3f & %20s & %20s & %10.2f & %10.2f \n'
                for iwvl,  amatch in enumerate(self.match):
                    chi = np.abs(self.Intensity[iwvl]/(2.*self.WghtFactor*amatch['predicted']))
                    pstring = pformat1L%(iwvl, amatch['obsWvl'], self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'], self.Intensity[iwvl]/amatch['predicted'])
                    outpt.write(pstring +'\n')
                    #
                    # now check line contributions
                    #
                    for jon,  anion in enumerate(amatch['ion']):
                        ionPrint = 0
                        for iline, awvl in enumerate(amatch['wvl'][jon]):
                            contrib = (amatch['intensity'][amatch['predictedLine'][jon][iline]]*self.Em).sum()/amatch['predicted']
                            if contrib > minContribution:
                                if ionPrint == 0:
                                    outpt.write(pformat2%(anion)+'\n')
                                    ionPrint = 1
                                outpt.write(pformat3L%(self.Wvl[iwvl], awvl, amatch['pretty1'][jon][iline], amatch['pretty2'][jon][iline],  self.Intensity[iwvl], amatch['predicted']))
                                outpt.write(dash +'\n')
                    outpt.write(dash +'\n')
        #
        # --------------------------------------------------------------------------
        #
    def predictPrint1d(self, minContribution=0.1, outfile=0, verbose=0):
        '''
        to predict the intensities of the observed lines from an emission measure
        the emission measure is already specified as self.Em which is an np array
        '''
        dash = ' -------------------------------------------------'
        pformat1 = ' %5i %7s %10.3f %10.2e %10.2e %10.3f %10.3f'
        pformat1s = ' %5s %7s %10s %10s %10s %10s %10s'
        pformat2 = '         %s'
        pformat3 = '        %10.3f %4i %4i %20s - %20s %5i %5i %7.3f'
        if self.Ndens > 1 and self.Ntemp == 1:
            idx = self.SearchData['best']['idx']
            dens = self.SearchData['best']['density']
            temp = self.SearchData['best']['temperature']
            emfit = self.SearchData['best']['emfit']
            em = self.SearchData['best']['em']
            print(' %5i %10.2e %10.2e %10.3f %10.2e'%(idx, dens, temp, emfit, em))
        else:
            print(' predictPrint1d does not know how to handle this')
            return
            # have just specified and Em by hand
        print(dash)
        pstring1 = pformat1s%('iwvl', 'IonS', 'wvl', 'Int',  'Pred', 'Int/Pred', 'chi')
        print(pstring1)
        print(dash)
        for iwvl,  amatch in enumerate(self.match):
            if amatch['predicted'] > 0. :
                chi = np.abs(self.Intensity[iwvl]/(2.*self.WghtFactor*amatch['predicted']))
                pstring = pformat1%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'], self.Intensity[iwvl]/amatch['predicted'], chi  )
            else:

                pstring = pformat1%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'], -1.)
            print(pstring)
            #
            # now check line contributions
            #
            for jon,  anion in enumerate(amatch['ion']):
                ionPrint = 0
                for iline, awvl in enumerate(amatch['wvl'][jon]):
                    contrib = (amatch['intensity'][amatch['predictedLine'][jon][iline]]*self.Em).sum()/amatch['predicted']
                    if contrib > minContribution:
                        if ionPrint == 0:
                            print(pformat2%(anion))
                            ionPrint = 1
                        print(pformat3%(awvl, amatch['lvl1'][jon][iline], amatch['lvl2'][jon][iline], amatch['pretty1'][jon][iline], amatch['pretty2'][jon][iline].ljust(20), amatch['lineIdx'][jon][iline], amatch['predictedLine'][jon][iline], contrib))
                        print(dash)
            print(dash)
        if outfile:
            with open(outfile, 'w') as outpt:
                try:
                    idx = self.SearchData['best']['idx']
                    dens = self.SearchData['best']['density']
                    temp = self.SearchData['best']['temperature']
                    emfit = self.SearchData['best']['emfit']
                    em = self.SearchData['best']['em']
                    print(' %5i %10.2e %10.2e %10.3f %10.2e'%(idx, dens, temp, emfit, em))
                except:
                    for emidx in self.EmIndices:
                        outpt.write(' %5i%12.2e %12.2e \n'%(emidx, self.Temperature[emidx], self.Em[emidx]))
                finally:
                    pass
                outpt.write(dash +'\n')
                for iwvl,  amatch in enumerate(self.match):
                    chi = np.abs(self.Intensity[iwvl]/(2.*self.WghtFactor*amatch['predicted']))
                    pstring = pformat1%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'], self.Intensity[iwvl]/amatch['predicted'], chi  )
                    outpt.write(pstring +'\n')
                    #
                    # now check line contributions
                    #
                    for jon,  anion in enumerate(amatch['ion']):
                        ionPrint = 0
                        for iline, awvl in enumerate(amatch['wvl'][jon]):
                            contrib = (amatch['intensity'][amatch['predictedLine'][jon][iline]]*self.Em).sum()/amatch['predicted']
                            if contrib > minContribution:
                                if ionPrint == 0:
                                    outpt.write(pformat2%(anion)+'\n')
                                    ionPrint = 1
                                outpt.write(pformat3%(awvl, amatch['lvl1'][jon][iline], amatch['lvl2'][jon][iline], amatch['pretty1'][jon][iline], amatch['pretty2'][jon][iline], amatch['lineIdx'][jon][iline], amatch['predictedLine'][jon][iline], contrib) +'\n')
                                outpt.write(dash +'\n')
                    outpt.write(dash +'\n')
        #
        # ------------------------------------------------------
        #
    def getChisq(self):
        '''
        return chisq
        '''
        weightedDiff,  msk = self.getWeightedDiff()
        chisq = np.sum(weightedDiff**2)
        return chisq, msk

    def getNormalizedChisq(self):
        '''
        return normalized chisq
        '''
        chisq,  mask = self.getChisq()
        normalChisq = chisq/float(self.Nobs - self.Nfree)
        return normalChisq
        #
        # -------------------------------------------------
        #
    def getWeightedDiff(self):
        '''
        to calculated the weighted difference of each of the intensities
        returns a 1D array
        '''
        self.predict()
        nwvl = len(self.match)
        weightedDiff = np.zeros(nwvl, 'float64')
        for iwvl, amatch in enumerate(self.match):
            if amatch['predicted'] > 0.:
                msk = False
                weightedDiff[iwvl] = (self.Intensity[iwvl] - amatch['predicted'])/(self.WghtFactor*self.Intensity[iwvl])
#                weightedDiff[iwvl] = (np.log10(self.Intensity[iwvl]) - np.log10(amatch['predicted']))/(self.WghtFactor*np.log10(amatch['predicted']))
            else:
                weightedDiff[iwvl] = 10.
                msk = True
        return weightedDiff, msk
        #
        #-----------------------------------------------------
        #
    def findMinMaxIndices(self, verbose=0):
        ''' to find the minimum and maximum indices where all match['intensitySum'] are
        greater than 0'''
        nT = len(self.Temperature)
        nlines = len(self.match)
        print(' n lines = %5i '%(nlines))
        minDx = 0
        maxDx = nT
        for iline, match in enumerate(self.match):
            if match['intensitySum'].sum() == 0.:
                print(' intensity Sum = 0 %5i  %10.3f'%(iline,  match['obsWvl']))
#            print('%5i %12.2e %12.2e'%(iline, match['intensitySum'].min(), match['intensitySum'].max()))
            else:
                idx1 = 0
                while match['intensitySum'][idx1] == 0.:
                    idx1 += 1
                idx2 = nT - 1
                while match['intensitySum'][idx2] == 0.:
                    idx2 -= 1
    #            print('min idx = %5i max idx =%5i'%(idx1,  idx2))
                if idx1 > minDx:
                    minDx = idx1
                if idx2 < maxDx:
                    maxDx = idx2
        print(' input # of temperatures %i'%(self.Temperature.size))
        print(' input temperature range %10.2e %10.2e '%(self.Temperature[0], self.Temperature[-1]))
        print(' self.MinIndex =     %5i      self.MaxIndex idx =%5i'%(minDx, maxDx))
        print(' min T =     %10.2e      max T = %10.2e'%(self.Temperature[minDx], self.Temperature[maxDx]))
        print(' log min T = %10.3f  log max T = %10.3f'%(np.log10(self.Temperature[minDx]), np.log10(self.Temperature[maxDx])))
        self.MinIndex = minDx
        self.MaxIndex = maxDx
        #
        # ------------------------------------------------------------------------
        #
    def fitFunc1t(self, em):
        '''
        the fitting function for the isothermal model to be called by leastsq

        Parameters
        ----------

        em:  number
            the log10 value of the emission measure

        Returns
        -------

        weighted chisquared:  float


        '''
        #
        self.emSet(em)
        self.predict()
        # change because getWeightedDiff also now returns a mask
        wd = self.getWeightedDiff()
        return wd[0]
        #
        # ------------------------------------------------------------------------
        #
    def fitFuncNt(self, value):
        '''
        the fitting function for the 1  (single temperature) temperature model to be called by leastsq
        '''
        self.emSet(value)
        self.predict()
        wd = self.getWeightedDiff()
        return wd[0]
        #
        # -------------------------------------------------------------
        #
    def fit1t(self, initialValue, maxfev=0):
        '''
        calls leastsq to fit the 1t (single temperature) model
        '''
        out = optimize.leastsq(self.fitFunc1t, np.asarray(initialValue, 'float64'), full_output=1, maxfev=maxfev)
        self.Leastsq = {'em':out[0], 'cov':out[1], 'info':out[2], 'message':out[3], 'ier':out[4]}
        #
        # -------------------------------------------------------------
        #
    def fitNt(self, initialValue, maxfev=0):
        '''
        calls leastsq to fit the 2d model
        '''
        out = optimize.leastsq(self.fitFuncNt, np.asarray(initialValue, 'float64'), full_output=1, maxfev=maxfev)
        self.Leastsq = {'em':out[0], 'cov':out[1], 'info':out[2], 'message':out[3], 'ier':out[4]}
        #
        # -----------------------------------------------------------
        #
    def search1dSpace(self, initialEm, indxlimits=None, verbose=0, log=0, maxfev=0):
        '''
        to conduct a brute force search over electron density for an isothermal-space and find the
        best fit to the em and density
        indxlimits give the range of indices to fit over
        can use self.MinIndex and self.MaxIndex+1
        initialEm = log value of the emission measure to begin the searching
        '''
        t1 = datetime.now()
        if self.Nobs <= 2.:
            print(' the number of observables %i is less than the number of parameters'%(self.Nobs))
        if indxlimits is None:
            if not hasattr(self, 'MinIndex'):
                self.findMinMaxIndices()
            indxlimits = [self.MinIndex, self.MaxIndex]
        self.Nfree = 2.

        emfit = []
        em = []
        idx = []
        chisq = []
        info = []
        searchDx = []
        maskedValues = []


        densSearched = np.asarray(self.EDensity[indxlimits[0]:indxlimits[1] + 1], 'float64')
        msk = np.ma.make_mask(np.ones((self.NTempDens)))
        chisq1d = np.ma.array(np.zeros((self.NTempDens), 'float64'), mask=msk)
        em1d = np.ma.array(np.zeros((self.NTempDens), 'float64'), mask=msk)
        if log:
            logname = 'search1d.raw'
            logfile = open(logname, 'w')
            logfile.write('%s \n'%(self.SpecData['filename']))
        logformat = '%5i %4i %10.2e %10.2e %10.2e  %10i\n'
        counter = 0
        pstring1 = '%5i %5i %10.2e %12.3e %10.2e %10.2e %10.2e '
        for idx1 in range(indxlimits[0], indxlimits[1] + 1, 1):
            searchDx1 = idx1 - indxlimits[0]
            searchDx.append(searchDx1)
            # kpd update
            self.emSetIndices(idx1)
            self.fit1t(initialEm, maxfev=maxfev)
            chisq1,  msk1 = self.getChisq()
            if not msk1:
                emfit.append(self.Leastsq['em'][0])
                chisq.append(chisq1)
                idx.append(idx1)
                info.append(self.Leastsq['info'])
                chisq1d.data[idx1] = chisq1
                chisq1d.mask[idx1] = False
                em1d.data[idx1] = self.Leastsq['em'][0]
                em1d.mask[idx1] = False
                if verbose:
                    em1 = 10.**(self.Leastsq['em'][0])
                    print(pstring1%(counter, idx1, self.EDensity[idx1], self.Leastsq['em'][0], em1, chisq[-1], self.Leastsq['info']['nfev']))
                elif log:
                    logfile.write(logformat%(counter, idx1, self.EDensity[idx1], self.Leastsq['em'][0], chisq[-1], self.Leastsq['info']['nfev']))
            else:
                maskedValues.append(idx1)
                if log:
                    logfile.write('msk at %5i %4i %10.2e \n'%(counter, idx1, chisq[-1]))
                # these values are already masked
                chisq1d.data[idx1] = self.getChisq()
                em1d.data[idx1] = self.Leastsq['em'][0]
            counter += 1
        if log:
            logfile.close()
    #
        em = 10.**np.asarray(emfit)
        self.SearchData = {'density':self.EDensity, 'densSearched':densSearched, 'emfit':emfit, 'em':em, 'idx':idx, 'chisq':chisq, 'minchisq':min(chisq), 'searchDx':searchDx, 'maskedValues':maskedValues, 'temperature':self.Temperature, 'message':'this contains all the data for the 1dSpace'}
        #
        print('min chisq = %10.3e '%(min(chisq)))
        # set things to best fit
        #
        gdx=[i for i,ch in enumerate(chisq) if ch == min(chisq)]
        #  kpd
        self.emSetIndices(idx[gdx[0]])
        self.emSet(emfit[gdx[0]])
        self.predict()
        em = 10.**emfit[gdx[0]]
        minChisq = min(chisq)
        reducedChisq = minChisq/float(self.Nobs - self.Nfree)
        #
        # key 'emfit' is the log value, 'em' is actual value
        self.SearchData['best'] = {'em':em, 'emfit':emfit[gdx[0]], 'chisq':chisq[gdx[0]], 'reducedChisq':reducedChisq, 'idx':idx[gdx[0]], 'density':self.EDensity[idx[gdx[0]]], 'temperature':self.Temperature[idx[gdx[0]]]}
        t2 = datetime.now()
        dt = t2 - t1
        print(' elapsed seconds = %12.3f'%(dt.seconds))
        #
        #-----------------------------------------------------
        #
    def search1tEmSpace(self, verbose=0):
        ''' to find the value of chisq as a function of Em with T = best-fit
        '''

        if not hasattr(self, 'SearchData'):
            print(' dem has not been searched yet')
            return
        if self.Nobs <= 2.:
            print(' the number of observables %10.2 is less than the number of parameters'%(self.Nobs, 2.*2.))
            return
        emChisq = []
        emfitSpace = self.SearchData['best']['emfit'] + np.linspace(-0.1,0.1,400) #*np.logspace(-3.,1.,21)   #n
        emSpace = 10.**(emfitSpace)
        self.emSetIndices(self.SearchData['best']['idx'])
        if verbose:
            print(' %5i  %10.2e'%(self.SearchData['best']['idx'], self.SearchData['best']['temperature']))
            print(' emfit %12.2e   em %12.2e '%(self.SearchData['best']['emfit'],self.SearchData['best']['em'] ))
        for idx,  anEm in enumerate(emfitSpace):
            #kpd
            self.emSet(anEm)
            chisq1,  msk1 = self.getChisq()
            if not msk1:
                emChisq.append(chisq1)
                em1 = np.log10(anEm)
                if verbose:
                    print(' %5i %12.3e %12.2e %12.2e '%(idx, anEm, em1, chisq1))
            else:
                    print('msk=1  %5i %12.3e %12.2e %12.2e '%(idx, anEm, em1, chisq1))

        #
        emChisq = np.asarray(emChisq, 'float64')
#        chisq1sig = emChisq.min()*np.ones_like(emChisq) + 2.3
#        good1sig = emChisq < emChisq.min() + 2.3
        print('min Em chisq = %10.3e'%(emChisq.min()))
#        em1sig = [emSpace[good1sig].min(), emSpace[good1sig].max()]
        #  will take care of confidence limits outside
        self.SearchData['EmError'] = {'em':emSpace, 'emfit':emfitSpace, 'chisq':emChisq, 'idx':self.SearchData['best']['idx'], 'temperature':self.Temperature[self.SearchData['best']['idx']]}

