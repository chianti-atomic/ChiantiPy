'''
classes and methods to analyze observations of astrophysical spectra
'''
import os
from datetime import date
from datetime import datetime
import pickle
import warnings

try:
    import multiprocessing as mp
except:
    print(' your version of Python does not support multiprocessing \n you will not be able to use mgofnt')
#
try:
    from ipyparallel import Client
except ImportError:
    warnings.warn("ipyparallel not found. You won't be able to use the ipymgofnt module")

import numpy as np
from scipy.interpolate import splrep, splev
import matplotlib.pyplot as plt
import ChiantiPy.core as ch
import ChiantiPy.tools.io as io
import ChiantiPy.tools.util as util
import ChiantiPy.tools.constants as const
import ChiantiPy.tools.data as chdata
import ChiantiPy.Gui as chGui
from ChiantiPy.base import ionTrails, specTrails    #
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
    # ---------------------------------------------------------
    #
def emPlot(matchDict, vs='T', loc='upper right', fs=10,  adjust=None, position='both', legend = True,  verbose=0):
    '''
    to plot line intensities divided by gofnt
    adjust is to provide an adjustment to the position of the labels
    position : one of 'both', 'right', 'left', or 'none'
    '''
    match = matchDict['match']
    temp = matchDict['Temperature']
    dens = matchDict['EDensity']
    nInt = len(match)
    if position == 'both':
        hpos = ['b']*nInt
    elif position == 'right':
        hpos = ['r']*nInt
    elif position == 'left':
        hpos = ['l']*nInt
    elif type(position) is list:
        if len(position) == nInt:
            hpos = position
            print(' position is the right length')
        else:
            print(' position is not the right length')
            return
    else:
        print(' position not understood')
        return

    if adjust is None:
        adjust = np.ones(nInt, np.float64)
#        print(' nInt = %5i'%(nInt))
    if not 'intensity' in match[0]:
        print(' must run mgofnt or gofnt first')
        return
    elif vs == 'T':
        if temp[0] == temp[-1]:
            ntemp = 1
        else:
            ntemp = temp.size
        if ntemp > 1:
            em = np.zeros((nInt, ntemp), 'float64')
            xvar = temp
        for idx in range(nInt):
            nonzed = match[idx]['intensitySum'] > 0.
            large = match[idx]['intensitySum'] >  match[idx]['intensitySum'].max()*0.01
            realgood = np.logical_and(nonzed, large)
            if realgood.sum() > 0 and match[idx]['obsIntensity'] > 0.:
                wvlstr = ' %6.2f'%(match[idx]['obsWvl'])
                em[idx][realgood] = match[idx]['obsIntensity']/match[idx]['intensitySum'][realgood]
                plt.loglog(xvar[realgood], em[idx][realgood], lw=2, label = wvlstr)
                if ntemp > 1:
#                    lblx = match[idx]['tmax']
                    lblx1 = temp[realgood][0]
                    lblx2 = temp[realgood][-1]
                else:
                    lblx2 = match[idx]['dmax']
                lbly1 = em[idx][realgood][0]
                lbly2 = em[idx][realgood][-1]
                plt.text(lblx1, lbly1, wvlstr.strip())
                plt.text(lblx2, lbly2, wvlstr.strip())
            else:
                print('no values for idx = %5i wvl = %8.2f'%(idx, match[idx]['wvl']))
                print(' nonzed = %i'%(nonzed.sum()))
                print('intensity = %10.2e'%(match[idx]['intensity']))
        if legend:
            plt.legend(loc=loc, fontsize=fs)
        plt.ylabel('Emission Measure (cm$^{-5}$)', fontsize=14)
        plt.xlabel('Temperature (K)',  fontsize=14)
    elif vs != 'T':
        if dens[0] == dens[-1]:
            ndens = 1
        else:
            ndens = len(dens)
        if ndens > 1:
            em = np.zeros((nInt, ndens), 'float64')
            xvar = dens
        for idx in range(nInt):
            nonzed = match[idx]['intensitySum'] > 0.
            large = match[idx]['intensitySum'] >  match[idx]['intensitySum'].max()*0.01
            realgood = np.logical_and(nonzed, large)
            if realgood.sum() > 0 and match[idx]['obsIntensity'] > 0.:
                wvlstr = ' %5.1f'%(match[idx]['obsWvl'])
                em[idx][realgood] = match[idx]['obsIntensity']/match[idx]['intensitySum'][realgood]
                plt.loglog(xvar[realgood], em[idx][realgood], lw=2, label=wvlstr)
                print(' %i  %5.1f'%(idx, match[idx]['obsWvl']))
                lblx = dens
                lbly = em[idx][realgood]*1.05
                if hpos[idx] =='b':
                    plt.text(lblx[0], lbly[0]*adjust[idx], wvlstr.strip())
                    plt.text(lblx[-1], lbly[-1]*adjust[idx], wvlstr.strip(), horizontalalignment='right')
                elif hpos[idx] == 'r':
                    plt.text(lblx[-1], lbly[-1]*adjust[idx], wvlstr.strip(), horizontalalignment='right')
                elif hpos[idx] == 'l':
                    plt.text(lblx[0], lbly[0]*adjust[idx], wvlstr.strip())

            else:
                print('no values for idx = %5i wvl = %8.2f'%(idx, match[idx]['wvl']))
                print(' nonzed = %i'%(nonzed.sum()))
                print('intensity = %10.2e'%(match[idx]['obsIntensity']))
        if legend:
            plt.legend(loc=loc, fontsize=fs)
        plt.ylabel('Emission Measure (cm$^{-5}$)', fontsize=14)
        plt.xlabel('Electron Density (cm$^{-3}$)',  fontsize=14)
    #
    # --------------------------------------------------------------------------
    #
def diffPrintMcF(specData,  matchDict, filename='diffPrintMcF.txt',  sort=None):
    '''
    calculates the weighted and straight differences between observed and predicted
    prints the values saves to a file from a PyMC run
    '''
    match = matchDict['match']
    temp = matchDict['Temperature']
    dens = matchDict['EDensity']
    if 'Emindices' in matchDict.keys():
        emIndices = matchDict['EmIndices']
    else:
        emIndices = None
        print('EmIndices missing from matchDict')
    if 'Em' in matchDict.keys():
        em = matchDict['Em']
    else:
        Em = None
        print('Em missing from matchDict')

    nMatch = len(match)
    if sort is None:
        sorter = range(nMatch)
    elif sort == 'ion':
        indexer = []
        for amatch in self.match:
            ionStr = amatch['exptIon']
            ionDict = util.convertName(ionStr)
            Z = ionDict['Z']
            stage = ionDict['Ion']
            number = Z*1000 + stage
            indexer.append(number)
        sorter = np.argsort(indexer)

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
        print(' %5s %12s %12s %12s'%('index',  'density',  'temperature',  'Em'''))
        outpt.write(' %5s %12s %12s %12s \n'%('index',  'density',  'temperature',  'Em'''))
        if emIndices is not None:
            for emidx in emIndices:
                print(' %5i %12.2e %12.2e %12.2e'%(emidx, dens[emidx],  temp[emidx], em[emidx]))
                outpt.write(' %5i %12.2e %12.2e %12.2e \n'%(emidx, dens[emidx],  temp[emidx], em[emidx]))
        print(dash)
        outpt.write(dash+'\n')
        print(' chi = abs(int/(wght*pred))  strDiff = (int - pred)/pred')
        outpt.write(' chi = abs(int/(2*wght*pred))  strDiff = (int - pred)/pred \n')
        print(sformat%('iwvl', 'ionS', 'wvl', 'intensity', 'predicted', 'int/pred', 'chi', 'relDev', 'stdDev'))
        outpt.write(fsformat%('iwvl', 'ionS', 'wvl', 'intensity', 'predicted', 'int/pred', 'chi', 'relDev',  'stdDev'))
        for iwvl in sorter:
            amatch = match[iwvl]
            if amatch['predicted'] > 0. :
                chi = np.abs(self.Intensity[iwvl]-amatch['predicted'])/(self.WghtFactor*self.Intensity[iwvl])
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
        normChisq = self.getNormalizedChisq()
        print('Normalized Chisq = %10.3f'%(normChisq))
        outpt.write('Normalized Chisq = %10.3f \n'%(normChisq))
        poor = diff > threeSig
        idx = np.arange(diff.size)
        pdx = idx[poor]
        for i in pdx:
            print('%5i %s %10.3f %10.3f'%(i, self.IonS[i],  self.Wvl[i], diff[i]))
            outpt.write('%5i %s %10.3f %10.3f \n'%(i, self.IonS[i],  self.Wvl[i], diff[i]))
    self.diffDict = {'diff':diff, 'wvl':self.Wvl, 'ionS':self.IonS, '3sig':threeSig, 'poor':poor}
#        self.SearchData['diff'] = diffDict
     #
    #-----------------------------------------------------
    #
def makeMatchPkl(specData, temp, dens, wghtFactor = 0.25,  abundanceName = None, minAbund=1.e-6, useMgofnt=1, verbose=0):
    '''
    input a data dictionary and instantiate a maker class,
    and run mgofnt and then make a pickle file
    to use multiprocessing, this needs to be run in an ipython console

    Parameters

    specData : dict
        the observed line intensities, wavelegths ...
    '''

    mydem = maker(specData, wghtFactor = wghtFactor, abundanceName =  abundanceName, minAbund=minAbund, verbose=verbose)
    mydem.makeMatch()
    if useMgofnt:
        mydem.mgofnt(temp, dens, verbose=1)
    else:
        mydem.gofnt(temp, dens)
#    mydem.predictor()
    filename = specData['filename']
    outname = os.path.splitext(filename)[0] + '-match-' + mydem.AbundanceName +'.pkl'
    mydem.PklName = outname
    print(' pickle name = %s'%(outname))
    mydem.dump(outname)
#    with open(outname, 'wb') as outpt:
#        pickle.dump(mydem.match, outpt)
    print(' should now move this to the working dir where it will be used')
    #
    #-----------------------------------------------------
    #
class maker(ionTrails,  specTrails):
    '''
    a class matching observed lines to lines in the CHIANTI database
    this class of MakerMc is for use with pymc3 MCMC modelling with
    continuous variables of temperature and emission measure
    '''
    def __init__(self, temperature, specData, elementList=[], ionList=[], allLines=False, abundanceName = None, minAbund=10., wghtFactor=None,  verbose=False):
        '''
        input a list of wavelengths and a wavelength difference
        find a list of predicted spectral lines for each wavelength
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
        elementList :  list
            a list of elements, such as fe, to be searched
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
        self.Temperature = temperature
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
        abundAll = chdata.Abundance[self.AbundanceName]['abundance']
        self.AbundAll = abundAll
        if wghtFactor is not None:
            self.WghtFactor = wghtFactor

        self.SpecData = specData
        self.ElementList = elementList
        self.IonList = ionList
        self.AllLines = allLines
        self.MinAbund = minAbund
        print(' minimum abundance = %10.2e'%(minAbund))
        self.MlInfo = io.masterListInfo()
        self.Intensity = specData['intensity']
        self.IonS = specData['ions']
        self.IonSet = set(specData['ions'])
        self.Dwvl = specData['dwvl']
        wvl = specData['wvl0']
        self.WvlRange = [min(wvl)-max(self.Dwvl), max(wvl)+max(self.Dwvl)]
        reduceNobs = 0.
        for anion in self.IonSet:
            same = specData['ions'].count(anion)
            reduceNobs += np.sqrt(float(same))
        self.ReduceNobs = reduceNobs
        self.Nobs = len(specData['wvl0'])
        print(' # of observables / reduce # = %i  %10.2f'%(self.Nobs,self.ReduceNobs))
        self.Wvl = specData['wvl0']
        #
#        self.WghtFactor = wghtFactor
#        print('wghtFactor = %8.3f'%(self.WghtFactor))
        #        self.SpecData['wghtFactor'] = []
#        if wghtFactor > 0.:
#            for iwvl in range(self.Nobs):
#                self.WghtFactor = self.SpecData['intStd'][iwvl]/self.SpecData['intensity'][iwvl] + wghtFactor
#        else:
#            for iwvl in range(self.Nobs):
#                self.WghtFactor = self.SpecData['intStd'][iwvl]/self.SpecData['intensity'][iwvl]
        #

    def makeMatch(self,  verbose=False):
        """ to match the CHIANTI lines with the input specdata
        uses ionTrails.ionGate to sort through ions
        """

        if 'intensity' in self.SpecData.keys():
            self.Intensity = self.SpecData['intensity']
#        masterlist = io.masterListRead()
#        matches = [{'ion':[], 'wvl':[]}]*len(wvl) - it won't work this way
        matches = []
        for iwvl in range(len(self.Wvl)):
            matches.append({'ion':[], 'wvl':[], 'lineIdx':[], 'wvldiff':[], 'lvl1':[], 'lvl2':[], 'pretty1':[], 'pretty2':[], 'predictedLine':[],'iPredictedLine':0, 'obsIntensity':self.Intensity[iwvl], 'exptIon':self.IonS[iwvl], 'obsWvl':self.Wvl[iwvl]})
        # use the ionList but make sure the ions are in the database

        self.ionGate(elementList = self.ElementList, ionList = self.IonList, minAbund=self.MinAbund, doLines=True, doContinuum=False, verbose = verbose)

        allLines = self.AllLines
        mlInfo = self.MlInfo
        for thision in self.Todo:
            if verbose:
#                print(' - - - - - - - - - - - ')
                print(' in makeMatchGate thision = %s'%(thision))
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

        if pDensity == 'default':
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
            if pDensity == 'default' and eDensity != None:
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
            if verbose:
                print(' using %s'%(someIon))
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
        self.ToProcess = []
        for someIon in ionList:
            # already know this ion is needed
            print(' someIon = %s'%(someIon))
            ionWorkerQ.put((someIon, temperature, density, self.AllLines))
            self.ToProcess.append(someIon)
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
    def tSet(self, value, verbose=0):
        '''
        to set the temperatures of the N temperature/density EM distribution
        '''
        if hasattr(self, 'Temperature'):
            self.Tguess = np.atleast_1d(value)
            self.Nfree = 2*len(value)
        else:
            print(' need to load pkl file or call gofnt/mgofnt first')
        #
        #------------------------------------------
        #
    def emSet(self, value):
        '''
        sets the EM values for a N temperature EM distribution
        '''
        emValue = np.atleast_1d(value)
        if hasattr(self, 'Tguess') and len(emValue) == len(self.Tguess):
            em = np.asarray(emValue, np.float64)
            self.Em = em
        else:
            print('in emSet, number of Tguess does not match that of emGuess ')
            if hasattr(self, 'Tguess '):
                print(' # of indices = %i, len(Tguess) = %i'%(len(value),len(self.Tguess)))
            else:
                print(' self.Tguess does not exist')
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
    def emMake(self, outName,  reference):
        """ to make an emission measure file
        outName does not need the suffix .em
        reference should be a list of references
        """
        if '.em' not in outName:
            outName += '.em'
#        print('writing file %s'%(outName))
        try:
            indices = self.EmIndices
        except:
            print(' the indices of for the prediction have not been set')
            return
#        for idx in indices:
#            print('T %10.3e  eD  %10.3e  EM  %10.3e'%(self.Temperature[idx], self.EDensity[idx], self.Em[idx]))

        with open(outName, 'w') as output:
            pformat = '%15.3e%15.3e%15.3e \n'
            for idx in indices:
                output.write(pformat%(self.Temperature[idx], self.EDensity[idx], self.Em[idx]))
            output.write(' -1\n')
            output.write('filename: %s\n'%(outName))
            for one in reference:
                output.write(one+'\n')
        #
        # ---------------------------------------------------------
        #
    def emPlot(self, vs='T', loc='upper right', fs=10,  adjust=None, position='both', label=True, legend = True, verbose=1):
        '''
        to plot line intensities divided by gofnt
        adjust is to provide an adjustment to the position of the labels
        position : one of 'both', 'right', 'left', or 'none'
        '''
        match = self.match
        temp = self.Temperature
        dens = self.EDensity
        nInt = len(match)
        if position == 'both':
            hpos = ['b']*nInt
        elif position == 'right':
            hpos = ['r']*nInt
        elif position == 'left':
            hpos = ['l']*nInt
        elif type(position) is list:
            if len(position) == nInt:
                hpos = position
                print(' position is the right length')
            else:
                print(' position is not the right length')
                return
        else:
            print(' position not understood')
            return

        if adjust is None:
            adjust = np.ones(nInt, np.float64)
    #        print(' nInt = %5i'%(nInt))
        if not 'intensity' in match[0]:
            print(' must run mgofnt or gofnt first')
            return
        elif vs == 'T':
            if temp[0] == temp[-1]:
                ntemp = 1
            else:
                ntemp = temp.size
#                ntemp = match[0]['intensity'].shape[0]
            if ntemp > 1:
                em = np.zeros((nInt, ntemp), 'float64')
                xvar = temp
            for idx in range(nInt):
                nonzed = match[idx]['intensitySum'] > 0.
                large = match[idx]['intensitySum'] >  match[idx]['intensitySum'].max()*0.01
                realgood = np.logical_and(nonzed, large)
                if realgood.sum() > 0 and match[idx]['obsIntensity'] > 0.:
                    wvlstr = ' %5.1f'%(match[idx]['obsWvl'])
                    em[idx][realgood] = match[idx]['obsIntensity']/match[idx]['intensitySum'][realgood]
                    if label:
                        plt.loglog(xvar[realgood], em[idx][realgood], lw=2, label = wvlstr)
                    else:
                        plt.loglog(xvar[realgood], em[idx][realgood], lw=2)
                    if ntemp > 1:
                        lblx1 = temp[realgood][0]
                        lblx2 = temp[realgood][-1]
                    else:
                        lblx1 = self.Temperature[realgood].min()
                        lbly1 = em[idx][0]
                    lbly1 = em[idx][realgood][0]
                    lbly2 = em[idx][realgood][-1]
                    plt.text(lblx1, lbly1, wvlstr.strip())
                    plt.text(lblx2, lbly2, wvlstr.strip())
                else:
                    print(' len of match[idx]wvl   %i'%(len(match[idx]['wvl'])))
                    if len(match[idx]['wvl']) != 0:
                        print('no values for idx = %5i wvl = %8.2f'%(idx, match[idx]['wvl']))
                        print(' nonzed = %i'%(nonzed.sum()))
                        print('intensity = %10.2e'%(match[idx]['intensity']))
            if legend:
                plt.legend(loc=loc, fontsize=fs)
            plt.ylabel('Emission Measure (cm$^{-5}$)', fontsize=14)
            plt.xlabel('Temperature (K)',  fontsize=14)
        elif vs != 'T':
            if dens[0] == dens[-1]:
                ndens = 1
            else:
                ndens = len(dens)
            if ndens > 1:
                em = np.zeros((nInt, ndens), 'float64')
                xvar = dens
            for idx in range(nInt):
                nonzed = match[idx]['intensitySum'] > 0.
                large = match[idx]['intensitySum'] >  match[idx]['intensitySum'].max()*0.01
                realgood = np.logical_and(nonzed, large)
                if realgood.sum() > 0 and match[idx]['obsIntensity'] > 0.:
                    wvlstr = ' %5.1f'%(match[idx]['obsWvl'])
                    em[idx][realgood] = match[idx]['obsIntensity']/match[idx]['intensitySum'][realgood]
                    if label:
                        plt.loglog(xvar[realgood], em[idx][realgood], lw=2, label=wvlstr)
                    else:
                        plt.loglog(xvar[realgood], em[idx][realgood], lw=2)
                    print(' %i  %5.1f'%(idx, match[idx]['obsWvl']))
                    lblx = dens
                    lbly = em[idx][realgood]*1.05
                    if hpos[idx] =='b':
                        plt.text(lblx[0], lbly[0]*adjust[idx], wvlstr.strip(), fontsize=14)
                        plt.text(lblx[-1], lbly[-1]*adjust[idx], wvlstr.strip(), horizontalalignment='right', fontsize=14)
                    elif hpos[idx] == 'r':
                        plt.text(lblx[-1], lbly[-1]*adjust[idx], wvlstr.strip(), horizontalalignment='right', fontsize=14)
                    elif hpos[idx] == 'l':
                        plt.text(lblx[0], lbly[0]*adjust[idx], wvlstr.strip(), fontsize=14)

                else:
                    print('no values for idx = %5i wvl = %8.2f'%(idx, match[idx]['wvl']))
                    print(' nonzed = %i'%(nonzed.sum()))
                    print('intensity = %10.2e'%(match[idx]['obsIntensity']))
            if legend:
                plt.legend(loc=loc, fontsize=fs)
            plt.ylabel('Emission Measure (cm$^{-5}$)', fontsize=14)
            plt.xlabel('Electron Density (cm$^{-3}$)',  fontsize=14)
            plt.tight_layout()
        #
        # --------------------------------------------------------------------------
        #
    def diffPrintMc(self, dir = '.', filename='diffPrintMc.txt',  sort=None):
        '''
        calculates the weighted and straight differences between observed and predicted
        prints the values saves to a file from a PyMC3 run
        also created a attribute self.Dict, a dict with the following keys:
        'wvl' = observed wavelength (A)
        'relDiff' = (I_obs - I_pred)/(I_obs)
        'ionS' the CHIANTI type name for an ion
        '''
        wghtFactor = self.WghtFactor
        today = date.today()
        thisday =today.strftime('%Y_%B_%d')
        nMatch = len(self.match)
        if sort is None:
            sorter = range(nMatch)
        elif sort == 'ion':
            indexer = []
            for amatch in self.match:
                ionStr = amatch['exptIon']
                ionDict = util.convertName(ionStr)
                Z = ionDict['Z']
                stage = ionDict['Ion']
                number = Z*1000 + stage
                indexer.append(number)
            sorter = np.argsort(indexer)

        wDiff = []
        intOverPred = []
        diffOverInt = []
        dash = ' -------------------------------------------------'
        pformat1 = ' %5i %7s %10.3f %10.2e %10.2e %10.3f %10.3f %10.3f %10.3f'
        pformat1a = ' %5i %7s %10.3f %10.2e %10.2e *****NID******'
        sformat  = ' %5s %7s %10s %10s %10s %10s %10s %10s %10s'
        fsformat  = ' %5s %7s %10s %10s %10s %10s %10s %10s %10s\n'
        with open(filename, 'w') as outpt:
            print(' today is %s'%(thisday))
            outpt.write(' today is %s \n'%(thisday))
            print(' WghtFactor = %10.3f'%(wghtFactor))
            outpt.write(' WghtFactor = %10.3f \n'%(wghtFactor))
            emIndices = self.EmIndices
            print(' %5s %12s %12s %12s'%('index',  'density',  'temperature',  'Em'''))
            outpt.write(' %5s %12s %12s %12s \n'%('index',  'density',  'temperature',  'Em'''))
            for emidx in emIndices:
                print(' %5i %12.2e %12.2e %12.2e'%(emidx, self.EDensity[emidx],  self.Temperature[emidx], self.Em[emidx]))
                outpt.write(' %5i %12.2e %12.2e %12.2e \n'%(emidx, self.EDensity[emidx],  self.Temperature[emidx], self.Em[emidx]))
            print(dash)
            outpt.write(dash+'\n')
            print(' chi = abs(int/(wght*pred))  strDiff = (int - pred)/pred')
            outpt.write(' chi = abs(int-pred)/(wght*pred))  strDiff = (int - pred)/pred \n')
            print(sformat%('', '', 'A', '', '', '', '', 'abs', 'abs' ))
            print(sformat%('iwvl', 'ionS', 'wvl', 'intensity', 'predicted', 'int/pred', 'chi', 'relDev', 'stdDev'))
            outpt.write(fsformat%('', '', 'A', '', '', '', '', 'abs', 'abs' ))
            outpt.write(fsformat%('iwvl', 'ionS', 'wvl', 'intensity', 'predicted', 'int/pred', 'chi', 'relDev',  'stdDev'))
            for iwvl in sorter:
                amatch = self.match[iwvl]
                if amatch['predicted'] > 0. :
                    chi = np.abs(self.Intensity[iwvl]-amatch['predicted'])/(wghtFactor*self.Intensity[iwvl])
                    wDiff.append(chi)
                    intOverPred.append(self.Intensity[iwvl]/amatch['predicted'] - 1.)
                    diffOverInt.append((self.Intensity[iwvl]-amatch['predicted'])/self.Intensity[iwvl])
                    pstring = pformat1%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'], self.Intensity[iwvl]/amatch['predicted'], chi, intOverPred[-1],  diffOverInt[-1] )
                else:
                    pstring = pformat1a%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'])
                print(pstring)
                outpt.write(pstring+'\n')
            intOverPredNp = np.asarray(intOverPred, np.float64)
            diffOverIntNp = np.asarray(diffOverInt, np.float64)
            print(' mean of relative Deviation = %10.3f 3*relDev = %10.3f stdDev = %10.3f\n'%(intOverPredNp.mean(), 3.*intOverPredNp.mean(), intOverPredNp.std()))
            outpt.write(' mean of rel Dev = %10.3f 3*rel Dev  = %10.3f stdDiff:  %10.3f \n'%(intOverPredNp.mean(), 3.*intOverPredNp.mean(), intOverPredNp.std()))
            print(' mean of intOverPred = abs(int/pred -1.) %10.3f'%(intOverPredNp.mean()))
            outpt.write(' mean of intOverPred = abs(int/pred -1.) %10.3f \n'%(intOverPredNp.mean()))
            threeSig = 3.*diffOverIntNp.std()
            print(' 3*std = %10.3f'%(threeSig))
            outpt.write(' 3*std = %10.3f \n'%(threeSig))
            chisq,  msk = self.getChisq()
            normChisq = self.getNormalizedChisq()
            print('           Chisq = %10.3f'%(chisq))
            print('Normalized Chisq = %10.3f'%(normChisq))
            outpt.write('           Chisq = %10.3f \n'%(chisq))
            outpt.write('Normalized Chisq = %10.3f \n'%(normChisq))
            poor = np.abs(diffOverIntNp) > threeSig
            idx = np.arange(nMatch)
            pdx = idx[poor]
            for i in pdx:
                print('%5i %s %10.3f %10.3f %10.3f'%(i, self.IonS[i],  self.Wvl[i], np.abs(intOverPredNp)[i],  diffOverIntNp[i]))
                outpt.write('%5i %s %10.3f %10.3f %10.3f\n'%(i, self.IonS[i],  self.Wvl[i], np.abs(intOverPredNp)[i], diffOverIntNp[i]))
        if sort is None:
            self.DiffMc = {'intOverPred':intOverPredNp, 'diffOverPred':diffOverIntNp, 'wvl':self.Wvl, 'ionS':self.IonS, '3sig':threeSig, 'poor':poor}
        else:
            self.DiffMcSort = {'intOverPred':intOverPredNp, 'diffOverPred':diffOverIntNp, 'wvl':self.Wvl, 'ionS':self.IonS, '3sig':threeSig, 'poor':poor}
        #
        # --------------------------------------------------------------------------
        #
    def predict(self):
        '''
        to predict the intensities of the observed lines from an emission measure
        the emission measure is already specified as self.Em which is an np array
        '''
        #

        nMatch = len(self.match)
        for iwvl, amatch in enumerate(self.match):


            bad = amatch['intensitySum'] <= 0.
            thismatch = np.ma.MaskedArray(amatch['intensitySum'],  dtype=np.float64, mask= bad)
            print('iwvl %i bad.sum %i'%(iwvl,  bad.sum()))


#            try:
#                tck = splrep(np.log(self.Temperature), np.log(thismatch),  s=0)
#                splpred  = splev(np.log(self.Tguess), tck,  der=0,  ext=1)
#                self.match[iwvl]cd['predicted'] = np.sum(splpred*self.Em)
#            except:
#                self.match[iwvl]['predicted'] = 0.
#                print(' error in predict, iwvl = %5i'%(iwvl))


            try:
                tck = splrep(self.Temperature, thismatch,  s=0)
                splpred  = splev(self.Tguess, tck,  der=0,  ext=1)
                self.match[iwvl]['predicted'] = np.sum(splpred*self.Em)
            except:
                self.match[iwvl]['predicted'] = 0.
                print(' error in predict, iwvl = %5i'%(iwvl))
        #
        # --------------------------------------------------------------------------
        #
    def predictPrint(self, minContribution=0.1, outfile='predictPrint.txt', sort=None, verbose=0):
        '''
        to predict the intensities of the observed lines from an emission measure
        the emission measure is already specified as self.Em which is an np array
        '''
        wghtFactor = self.WghtFactor
        nMatch = len(self.match)
        if sort is None:
            sorter = range(nMatch)
        elif sort == 'ion':
            indexer = []
            for amatch in self.match:
                ionStr = amatch['exptIon']
                ionDict = util.convertName(ionStr)
                Z = ionDict['Z']
                stage = ionDict['Ion']
                number = Z*1000 + stage
                indexer.append(number)
            sorter = np.argsort(indexer)

        dash = ' -------------------------------------------------'
        pformat1 = ' %5i %7s %10.3f %10.2e %10.2e %10.3f %10.3f'
        pformat1s = ' %5s %7s %10s %10s %10s %10s %10s'
        pformat2 = '         %s'
        pformat3 = '        %10.3f %4i %4i %20s - %20s %5i %5i %7.3f'
        pformat3s = '        %10s %4s %4s %20s - %20s %5s %5s %s'
        if outfile:
            with open(outfile, 'w') as outpt:
#                try:
#                    emIndices = self.EmIndices
#                    for emidx in emIndices:
#                        print(' %5i %12.2e %12.2e %12.2e'%(emidx, self.EDensity[emidx],  self.Temperature[emidx], self.Em[emidx]))
#                        outpt.write(' %5i %12.2e %12.2e %12.2e\n'%(emidx, self.EDensity[emidx],  self.Temperature[emidx], self.Em[emidx]))
#                except:
#                    print(' attribute EmIndices has not been created')
#                    outpt.write(' attribute EmIndices has not been created \n')
#                    return
                print('matchPkl %s '%(self.MatchName))
                outpt.write('matchPkl %s \n'%(self.MatchName))
                print('wghtFactor %10.3f'%(wghtFactor))
                outpt.write('wghtFactor %10.3f \n'%(wghtFactor))
                print(dash)
                outpt.write(dash + '\n')
                pstring1 = pformat1s%('iwvl', 'IonS', 'wvl', 'Int',  'Pred', 'Int/Pred', 'chi')
                print(pstring1)
                outpt.write(pstring1 +'\n')
                pstring3 = pformat3s%('wvl', 'lvl1', 'lvl2', 'lower', 'upper', 'lineIdx',  'predLine', 'contribution')
                print(pstring3)
                outpt.write(pstring3+'\n')
                print(dash)
                outpt.write(dash +'\n')
                for iwvl in sorter:
                    amatch = self.match[iwvl]
                    if amatch['predicted'] > 0. :
                        chi = (np.abs(self.Intensity[iwvl]-amatch['predicted'])/(wghtFactor*self.Intensity[iwvl]))
                        pstring = pformat1%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'], self.Intensity[iwvl]/amatch['predicted'], chi  )
                    else:

                        pstring = pformat1%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'], -1.)
                    print(pstring)
                    outpt.write(pstring +'\n')
                    #
                    # now check line contributions
                    #
                    for jon,  anion in enumerate(amatch['ion']):
                        if verbose:
                            print(' ion %s  jon %i'%(anion,  jon))
                        ionPrint = 0
                        for iline, awvl in enumerate(amatch['wvl'][jon]):
                            thisIntensity1 = amatch['intensity'][amatch['predictedLine'][jon][iline]]
                            bad = thisIntensity1 <= 0.
                            thismatch = np.ma.MaskedArray(thisIntensity1,  dtype=np.float64, mask= bad)
                            tck = splrep(self.Temperature, thismatch,  s=0)
                            splpred  = splev(self.Tguess, tck,  der=0,  ext=1)

#                            contrib = (amatch['intensity'][amatch['predictedLine'][jon][iline]]*self.Em).sum()/amatch['predicted']
                            contrib = np.sum(splpred*self.Em)
                            if contrib > minContribution:
                                if ionPrint == 0:
                                    print(pformat2%(anion))
                                    outpt.write(pformat2%(anion)+'\n')
                                    ionPrint = 1
                                print(pformat3%(awvl, amatch['lvl1'][jon][iline], amatch['lvl2'][jon][iline], amatch['pretty1'][jon][iline], amatch['pretty2'][jon][iline].ljust(20), amatch['lineIdx'][jon][iline], amatch['predictedLine'][jon][iline], contrib))
                                outpt.write(pformat3%(awvl, amatch['lvl1'][jon][iline], amatch['lvl2'][jon][iline], amatch['pretty1'][jon][iline], amatch['pretty2'][jon][iline].ljust(20), amatch['lineIdx'][jon][iline], amatch['predictedLine'][jon][iline], contrib)+'\n')
                                print(dash)
                                outpt.write(dash +'\n')
                    print(dash)
                outpt.write(dash +'\n')
                print('matchPkl %s \n'%(self.MatchName))
                print('wghtFactor %10.3f \n'%(wghtFactor))
                pstring1 = pformat1s%('iwvl', 'IonS', 'wvl', 'Int',  'Pred', 'Int/Pred', 'chi')
                outpt.write(pstring1 +'\n')
                outpt.write(pstring3 + '\n')
                for iwvl,  amatch in enumerate(self.match):
                    chi = (np.abs(self.Intensity[iwvl]-amatch['predicted'])/(2.*wghtFactor*self.Intensity[iwvl]))
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
                    chi = (np.abs(self.Intensity[iwvl]-amatch['predicted'])/(2.*wghtFactor*self.Intensity[iwvl]))
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
    def loadMatch(self, filename):
        """ to open a pickle file, return the match data and make it an attribute
        """
        with open(filename, 'rb') as inpt:
            matchDict = pickle.load(inpt)
        self.match = matchDict['match']
        self.Temperature = matchDict['Temperature']
        self.EDensity = matchDict['EDensity']
        self.Density = matchDict['EDensity']
        self.Ndens = matchDict['Ndens']
        self.Ntemp = matchDict['Ntemp']
        self.NTempDens = matchDict['NTempDens']
        if 'EmIndices' in matchDict.keys():
            self.EmIndices = matchDict['EmIndices']
        else:
            print(' EmIndices not in matchDict')
        if 'Em' in matchDict.keys():
            self.Em = matchDict['Em']
        else:
            print('Em not in matchDict')
        if 'XUVTOP' in matchDict.keys():
            self.XUVTOP = matchDict['XUVTOP']
        if 'chiantiVersion' in matchDict.keys():
            self.ChiantiVersion = matchDict['chiantiVersion']
        self.MatchName = filename
        #
        #-----------------------------------------------------
        #
    def dumpMatch(self, filename):
        """to save the attribute match to a pickle file
        """
        matchDict={'match':self.match, 'Temperature':self.Temperature, 'EDensity':self.EDensity, 'Ndens':self.Ndens,
            'Ntemp':self.Ntemp, 'NTempDens':self.NTempDens, 'MinAbund':self.MinAbund}
        if hasattr(self, 'EmIndices'):
            matchDict['EmIndices'] = self.EmIndices
        if hasattr(self, 'Em'):
            matchDict['Em'] = self.Em
        if hasattr(self, 'Nfree'):
            matchDict['Nfree'] = self.Nfree
        if hasattr(self, 'WghtFactor'):
            matchDict['WghtFactor'] = self.WghtFactor
        else:
            print('wghtfactor not available')
        if 'XUVTOP' in self.SpecData.keys():
            matchDict['XUVTOP'] = self.SpecData['XUVTOP']
        if 'chiantiVersion' in self.SpecData.keys():
            matchDict['chiantiVersion'] = self.SpecData['chiantiVersion']
        with open(filename, 'wb') as outpt:
            pickle.dump(matchDict, outpt)
