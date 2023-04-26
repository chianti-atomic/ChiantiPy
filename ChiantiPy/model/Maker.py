'''
classes and methods to analyze observations of astrophysical spectra
'''
import os
from datetime import date
from datetime import datetime
import pickle
import multiprocessing as mp


import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import ChiantiPy.core as ch
import ChiantiPy.tools.io as io
import ChiantiPy.tools.util as util
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


    Keyword Arguments
    -----------------

    vs: `str`, either 'T', or 'D'
        whether to plot the emission measure vs temperature or density

    loc:  `str`
        matplotlib argument for plt.legend

    fs:  `int`
        the fontsize for the legend

    adjust:  `list`
        a list of multiplicative adjustments to the labels to the plot lines
        must be the same length as the number of lines

    position:  `str`
        where the labels to the lines should be placed, `both` for both ends, `left` for the left size only, 'right' for the right side only, or None for no labels

    label:  `bool`
        whether to apply

    legend:  `bool`
        whether to include a matplotlib legend

    fontsize:  `int`
        fontsize for the matplotlib xlabel and ylabel

    tscale:  `float`
        scale the temperature by dividing by tscale

    verbose : `bool`
        if True, additional output is sent to the terminal


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
    elif position == 'none':
        hpos = ['n']*nInt
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
                if hpos[idx] == 'l':
                    plt.text(lblx1, lbly1, wvlstr.strip())
                elif hpos[idx] == 'r':
                    plt.text(lblx2, lbly2, wvlstr.strip())
                elif hpos[idx] == 'b':
                    plt.text(lblx1, lbly1, wvlstr.strip())
                    plt.text(lblx2, lbly2, wvlstr.strip())
                elif hpos[idx] == 'n':
                    pass
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
    #-----------------------------------------------------
    #
class maker(ionTrails,  specTrails):
    '''
    a class matching observed lines to lines in the CHIANTI database


    Parameters
    ----------
    specData : dict
        contains the following keys
        intensity - a list of observed line intensities
        wvlObs - a list of observed wavelengths, usually in Angstroms
        dwvl the expected wavelength different between the observed wvl and CHIANTI

    Keyword Arguments
    -----------------

    temperature: `float`, `list`, `ndarray`
        the temperature(s) in K

    eDensity: float, ndarray
        eDensity: electron density in :math:`\mathrm{cm^{-3}}`

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
    def __init__(self, specData, temperature=None, eDensity=None, elementList=[], ionList=[],
        allLines=False, abundanceName = None, minAbund=10., wghtFactor=None,  verbose=False):
        '''
        input a list of wavelengths and a wavelength difference
        find a list of predicted spectral lines for each wavelength
        all = 0 -> only previously observed wavelengths are used
        dwvl = wavelength difference between CHIANTI and the observed lines
            if == 0 then use use that in the observed data

        '''
        #
        # --------------------------------------------------------------------------------
        #
        if temperature is not None:
            self.Temperature = np.array(temperature, np.float64)
        if eDensity is not None:
            self.EDensity = np.array(eDensity, np.float64)
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
        self.Abundance = abundAll
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
        self.Nions = len(self.IonSet)
        self.Dwvl = specData['dwvl']
        wvl = specData['wvl0']
        try:
            self.WvlRange = [min(wvl)-max(self.Dwvl), max(wvl)+max(self.Dwvl)]
        except:
            if not hasattr(self, 'WvlRange'):
                print(' need to set WvlRange attribute')
                return
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


        Keyword Arguments
        -----------------

        verbose : `bool`
            if True, additional output is sent to the terminal

        """

        if 'intensity' in self.SpecData.keys():
            self.Intensity = self.SpecData['intensity']
#        masterlist = io.masterListRead()
#        matches = [{'ion':[], 'wvl':[]}]*len(wvl) - it won't work this way
        matches = []
        for iwvl in range(len(self.Wvl)):
            matches.append({'ion':[], 'wvl':[], 'lineIdx':[], 'wvldiff':[], 'lvl1':[], 'lvl2':[],
                'pretty1':[], 'pretty2':[], 'predictedLine':[],'iPredictedLine':0,
                'obsIntensity':self.Intensity[iwvl], 'exptIon':self.IonS[iwvl], 'obsWvl':self.Wvl[iwvl]})
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
        self.Match = matches


    def matchPrint(self, filename='matchPrint.txt', sort=None, verbose=False):
        '''
        to print out the data for the matches to the observed lines
        sort can be 'wvl' or 'ion', otherwise, there is no sorting done
        does not require predicted values


        Keyword Arguments
        -----------------


        filename:  `str`
            the filename where the text should be output

        sort:  `str`, can be `wvl`, `ion`, or None
            whether the output should be sorted by `wvl` or `ion` or not

        verbose : `bool`
            if True, additional output is sent to the terminal

        '''

        minContribution = 0.1

        cwd = os.getcwd()
        wghtFactor = 0.
        nMatch = len(self.Match)
        if verbose:
            print('nMatch:  %i'%(nMatch))
        if sort is None:
            sorter = range(nMatch)
        elif sort == 'ion':
            indexer = []
            for amatch in self.Match:
                ionStr = amatch['exptIon']
                ionDict = util.convertName(ionStr)
                Z = ionDict['Z']
                stage = ionDict['Ion']
                number = Z*1000 + stage
                indexer.append(number)
            sorter = np.argsort(indexer)
        elif sort == 'wvl':
            indexer = []
            for amatch in self.Match:
                wvlObs = amatch['obsWvl']
                indexer.append(wvlObs)
            sorter = np.argsort(indexer)

        dash = ' -------------------------------------------------'

        pformat1 = ' %5i %7s %10.3f %10.2e'
        pformat1s = ' %5s %7s %10s %10s'
        pformat2 = '         %s'
        pformat3 = '        %10.3f %10.3f %10.3f %4i %4i %20s - %20s %5i %5i'
        pformat3s = '       %10s %10s %10s %4s %4s %20s - %20s %5s %5s'
        if filename:
            with open(filename, 'w') as outpt:
                print(' cwd:  %s'%(cwd))
                outpt.write(' cwd:  %s \n'%(cwd))
                if hasattr(self, 'MatchName'):
                    print('matchPkl %s '%(self.MatchName))
                    outpt.write('matchPkl %s \n'%(self.MatchName))
                print('wghtFactor %10.3f'%(wghtFactor))
                outpt.write('wghtFactor %10.3f \n'%(wghtFactor))
                print(dash)
                outpt.write(dash + '\n')
                pstring1 = pformat1s%('iwvl', 'IonS', 'wvl', 'Int')
                print(pstring1)
                outpt.write(pstring1 +'\n')
                pstring3 = pformat3s%('wvl_obs','wvl', 'diff', 'lvl1', 'lvl2', 'lower',
                    'upper', 'lineIdx',  'predLine')
                print(pstring3)
                outpt.write(pstring3+'\n')
                print(dash)
                outpt.write(dash +'\n')
                for iwvl in sorter:
                    amatch = self.Match[iwvl]
                    pstring = pformat1%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl])
                    print(pstring)
                    outpt.write(pstring +'\n')
                    #
                    # now check line contributions
                    #
                    for jon,  anion in enumerate(amatch['ion']):
                        ionPrint = 0
                        for iline, awvl in enumerate(amatch['wvl'][jon]):
                            contrib = 1.
                            if contrib > minContribution:
                                if ionPrint == 0:
                                    print(pformat2%(anion))
                                    outpt.write(pformat2%(anion)+'\n')
                                    ionPrint = 1
                                dwvl = awvl - self.Wvl[iwvl]
                                print(pformat3%(self.Wvl[iwvl], awvl, dwvl,  amatch['lvl1'][jon][iline],
                                    amatch['lvl2'][jon][iline], amatch['pretty1'][jon][iline],
                                    amatch['pretty2'][jon][iline].ljust(20),
                                    amatch['lineIdx'][jon][iline], amatch['predictedLine'][jon][iline]))
                                outpt.write(pformat3%(self.Wvl[iwvl], awvl, dwvl, amatch['lvl1'][jon][iline],
                                    amatch['lvl2'][jon][iline], amatch['pretty1'][jon][iline],
                                    amatch['pretty2'][jon][iline].ljust(20), amatch['lineIdx'][jon][iline],
                                    amatch['predictedLine'][jon][iline])+'\n')
                                print(dash)
                                outpt.write(dash +'\n')
                    print(dash)
                outpt.write(dash +'\n')

    def argCheck(self, temperature=None, eDensity=None, pDensity='default', verbose=False):
        ''' to check the compatibility of the three arguments
        and put them into numpy arrays of atleast_1d

        Keyword Arguments
        -----------------

        temperature: `float`, `list`, `ndarray`
            the temperature(s) in K

        eDensity: float, ndarray
            eDensity: electron density in :math:`\mathrm{cm^{-3}}`

        pDensity: `str`, `float`, `ndarray`
            pDensity: proton density in :math:`\mathrm{cm^{-3}}` defaults to 'default'

        verbose : `bool`
            if True, additional output is sent to the terminal

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

        Parameters
        -----------------

        temperature: `float`, `list`, `ndarray`
            the temperature(s) in K

        density: `float`, `list`, `ndarray`
            density: electron density in :math:`\mathrm{cm^{-3}}`


        '''
        t1 = datetime.now()
        self.XUVTOP = os.environ['XUVTOP']
        #  quick and dirty
        #   em = None,
        self.argCheck(temperature=temperature, eDensity=density, pDensity=None,  verbose=0)
        temperature = self.Temperature
        density = self.EDensity

        nTempDens = self.NTempDens
        if verbose:
            print(' temperature size:  %5i'%(self.Temperature.size))
            print(' density     size:  %5i'%(self.EDensity.size))
        ionList = []
        for iwvl, amatch in enumerate(self.Match):
            for someIon in amatch['ion']:
#                print 'someIon = ', someIon
                if someIon not in ionList:
                    ionList.append(someIon)
        for iwvl, amatch in enumerate(self.Match):
            nPredictedLines = 0
            for awvl in amatch['wvl']:
                nPredictedLines +=  len(awvl)
            self.Match[iwvl]['intensity'] = np.zeros((nPredictedLines,  nTempDens), 'float64')

        self.ionList = ionList
        for iwvl in range(len(self.Match)):
            self.Match[iwvl]['intensitySum'] = np.zeros(nTempDens, 'float64')
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
                for iwvl, amatch in enumerate(self.Match):
    #                self.Match[iwvl]['intensitySum'] = np.zeros(nTempDens, 'float64')
                    #  this is data for each line
                    if someIon in amatch['ion']:
                        kon = amatch['ion'].index(someIon)
    #                for kon, anIon in enumerate(amatch['ion']):
    #
        #                print ' linIdx = ', amatch['lineIdx'], amatch
                        predictedLine = []
                        for aline in amatch['lineIdx'][kon]:
                            #print(' %s   %6i  %12.3f '%( someIon, aline, thisIon.Intensity['wvl'][aline]))
                            self.Match[iwvl]['intensitySum'] += intensity[:, aline]
                            iPredictedLine = self.Match[iwvl]['iPredictedLine']
                            predictedLine.append(iPredictedLine)
                            self.Match[iwvl]['intensity'][iPredictedLine] = intensity[:, aline]
                            self.Match[iwvl]['iPredictedLine'] += 1
                        self.Match[iwvl]['predictedLine'][kon] = predictedLine
                self.Tmax = np.zeros_like(self.Wvl)
#                self.Dmax = np.zeros_like(self.Wvl)
                for iwvl, amatch in enumerate(self.Match):
                    gfun = amatch['intensitySum']
                    peak = gfun == gfun.max()
#                    if len(temperature) > 1:
                    self.Match[iwvl]['tmax'] = temperature[peak]
#                    elif len(density) >1:
#                        self.Match[iwvl]['dmax'] = density[peak]
            else:
               print(' IoneqOne not available for %s'%(someIon))
        t2 = datetime.now()
        dt = t2 - t1
        print(' elapsed seconds = %12.3f'%(dt.seconds))
        #for iwvl, amatch in enumerate(self.Match):
            #gfun = amatch['intensitySum']
            #peak = gfun == gfun.max()
            #self.Match[iwvl]['tmax'] = temperature[peak]
        #
        # ---------------------------------------------------------------------
        #
    def mgofnt(self, temperature, density, proc=6,  timeout=0.1, verbose=False):
        '''
        calculate the gofnt function for each of the matched lines
        this is the multiprocessing version
        does each ion only once

        Parameters
        -----------------

        temperature: `float`, `list`, `ndarray`
            the temperature(s) in K

        density: `float`, `list`, `ndarray`
            density: electron density in :math:`\mathrm{cm^{-3}}`

        Keyword Arguments
        -----------------

        proc:  `int`
            the number of cores to be used

        timeout:  'float'
            may not actually be necessary

        verbose : `bool`
            if True, additional output is sent to the terminal


        '''
        t1 = datetime.now()
        self.XUVTOP = os.environ['XUVTOP']
        #  quick and dirty
        self.argCheck(temperature=temperature, eDensity=density, pDensity=None,  verbose=0)

        temperature = self.Temperature
        density = self.EDensity
        nTempDens = self.NTempDens

        for iwvl, amatch in enumerate(self.Match):
            nPredictedLines = 0
            for awvl in amatch['wvl']:
                nPredictedLines +=  len(awvl)
            self.Match[iwvl]['intensity'] = np.zeros((nPredictedLines,  nTempDens), 'float64')

        ionList = []
        for iwvl, amatch in enumerate(self.Match):
            for someIon in amatch['ion']:
#                print 'someIon = ', someIon
                if someIon not in ionList:
                    ionList.append(someIon)

        self.ionList = ionList

        for iwvl in range(len(self.Match)):
            self.Match[iwvl]['intensitySum'] = np.zeros(nTempDens, 'float64')
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

            for iwvl, amatch in enumerate(self.Match):
#                self.Match[iwvl]['intensitySum'] = np.zeros(nTempDens, 'float64')
                #  this is data for each line
                if someIon in amatch['ion']:
                    kon = amatch['ion'].index(someIon)
#                    if verbose:
#                        print('kon: %5i using %s'%(kon, someIon))
#                        if 'errorMessage' in out[1].keys():
#                            print(' in mgofnt %s'%(out[1]['errorMessage']))

                    predictedLine = []
                    for aline in amatch['lineIdx'][kon]:
#                        if verbose:
#                            print(' ion, lineIdx = %s %5i'%(someIon, aline))
                        self.Match[iwvl]['intensitySum'] += intensity[:, aline]
                        iPredictedLine = self.Match[iwvl]['iPredictedLine']
                        predictedLine.append(iPredictedLine)
                        self.Match[iwvl]['intensity'][iPredictedLine] = intensity[:, aline]
                        self.Match[iwvl]['iPredictedLine'] += 1
                    self.Match[iwvl]['predictedLine'][kon]=predictedLine
        #
        self.Tmax = np.zeros_like(self.Wvl)
        nT = self.Temperature.size
        trange = np.arange(nT)
        minIdx = []
        maxIdx = []
        for iwvl, amatch in enumerate(self.Match):
            gfun = amatch['intensitySum']
            peak = gfun == gfun.max()
            if len(temperature) > 1:
                self.Match[iwvl]['tmax'] = temperature[peak]
            elif len(density) >1:
                self.Match[iwvl]['dmax'] = density[peak]
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
    def emSetIndices(self, indices, add=0.,  verbose=0):
        '''
        to set the indices of the N temperature/density EM distribution
        can increase the number of paramaters if additional parameters have been used


        Parameters
        ----------

        indices:  `list`, `ndarray`
            the indices of the temperature/density arrays for which a set of intensities will be predicted


        Keyword Arguments
        -----------------

        add:  `float`
            to increase the number of parameters used in the calculation of the reduced chi-squared

        verbose : `bool`
            if True, additional output is sent to the terminal

        '''
        if hasattr(self, 'Temperature'):
            self.EmIndices = np.atleast_1d(indices)
            self.Nparams = 2.*self.EmIndices.size + add
            self.NT = self.EmIndices.size
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


        Parameters
        ----------

        value:  `list`, `ndarray`
            the values of the emission measure to be used when the intensities are predicted

        '''
        emValue = np.atleast_1d(value)
        if hasattr(self, 'EmIndices') and emValue.size == self.EmIndices.size:
            em = np.zeros_like(self.Temperature)
            emLog = np.zeros_like(self.Temperature)
            for i, adx in enumerate(self.EmIndices):
                em[adx] = 10.**emValue[i]
                emLog[adx] = emValue[i]
            self.Em = em
            self.EmLog = emLog
        else:
            print('in emNtexp, either EmIndices not set or no. of values != no. of indices ')
            if hasattr(self, 'EmIndices'):
                print(' # of indices = %i, EmIndices.size = %i'%(value.size, self.EmIndices.size))
            else:
                print(' self.EmIndices does not exist')
                print('# of indices = %i'%(emValue.size))
        #
        # ---------------------------------------------------------
        #
    def emFitPlot(self):
        '''
        to plot the emission measures derived from a chi-squared search over temperature
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
            plt.loglog([float(self.Temperature.min()), float(self.Temperature.max())], [em, em],
                '-k', linewidth=2 )
        #
        # ---------------------------------------------------------
        #
    def emWrite(self, filename, directory, reference):
        """ to make a CHIANTI style emission measure file
        outName does not need the suffix .em
        reference should be a list of references


        Parameters
        ----------

        filename:  `str`
            the name of the em file to be produced

        reference:  `list`
            a list of strings providing a reference at the end of the em file


        """
        if '.em' not in filename:
            filename += '.em'
#        print('writing file %s'%(outName))
        fullFileName = os.path.join(directory,  filename)
        try:
            indices = self.EmIndices
        except:
            print(' the indices of for the prediction have not been set')
            return
#        for idx in indices:
#            print('T %10.3e  eD  %10.3e  EM  %10.3e'%(self.Temperature[idx], self.EDensity[idx], self.Em[idx]))

        with open(fullFileName, 'w') as output:
            pformat = '%15.3e%15.3e%15.3e \n'
            for idx in indices:
                output.write(pformat%(self.Temperature[idx], self.EDensity[idx], self.Em[idx]))
            output.write(' -1\n')
            output.write('filename: %s\n'%(filename))
            for one in reference:
                output.write(one+'\n')
        #
        # ---------------------------------------------------------
        #
    def emPlot(self, vs='T', loc='upper right', fs=10,  adjust=None, position='both', label = True, legend = True, fontsize=16, tscale=1.,   verbose=True):
        '''
        to plot line intensities divided by gofnt
        adjust is to provide an adjustment to the position of the labels
        position : one of 'both', 'right', 'left', or 'none'

        Keyword Arguments
        -----------------

        vs: `str`, either 'T', or 'D'
            whether to plot the emission measure vs temperature or density

        loc:  `str`
            matplotlib argument for plt.legend

        fs:  `int`
            the fontsize for the legend

        adjust:  `list`
            a list of multiplicative adjustments to the labels to the plot lines
            must be the same length as the number of lines

        position:  `str`
            where the labels to the lines should be placed, `both` for both ends, `left` for the left size only, 'right' for the right side only, or None for no labels

        label:  `bool`
            whether to apply

        legend:  `bool`
            whether to include a matplotlib legend

        fontsize:  `int`
            fontsize for the matplotlib xlabel and ylabel

        tscale:  `float`
            scale the temperature by dividing by tscale

        verbose : `bool`
            if True, additional output is sent to the terminal

       '''
        match = self.Match
        temp = self.Temperature
        dens = self.EDensity
        nInt = len(match)
        if position == 'both':
            hpos = ['b']*nInt
        elif position == 'right':
            hpos = ['r']*nInt
        elif position == 'left':
            hpos = ['l']*nInt
        elif position is None:
            hpos = ['n']*nInt
        elif type(position) is list:
            if len(position) == nInt:
                hpos = position
                print(' position is the right length')
            else:
                print(' position is not the right length')
                return
        else:
            print(' position not understood: should be left, right, both, None')
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
                xvar = temp/tscale
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
                    if hpos[idx] == 'l':
                        plt.text(lblx1, lbly1, wvlstr.strip())
                    elif hpos[idx] == 'r':
                        plt.text(lblx2, lbly2, wvlstr.strip())
                    elif hpos[idx] == 'b':
                        plt.text(lblx2, lbly2, wvlstr.strip())
                        plt.text(lblx1, lbly1, wvlstr.strip())
                else:
                    print(' len of match[idx]wvl   %i'%(len(match[idx]['wvl'])))
                    if len(match[idx]['wvl']) != 0:
                        print('no values for idx = %5i wvl = %8.2f'%(idx, match[idx]['wvl']))
                        print(' nonzed = %i'%(nonzed.sum()))
                        print('intensity = %10.2e'%(match[idx]['intensity']))
            if legend:
                plt.legend(loc=loc, fontsize=fontsize)
            plt.ylabel('Emission Measure (cm$^{-5}$)', fontsize=fontsize)
            plt.xlabel('Temperature (K)',  fontsize=fontsize)
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
            plt.ylabel('Emission Measure (cm$^{-5}$)', fontsize=fontsize)
            plt.xlabel('Electron Density (cm$^{-3}$)',  fontsize=fontsize)
            plt.tight_layout()
        #
        # ---------------------------------------------------------
        #
    def emPlotObj(self, vs='T', loc='upper right', fs=10,  adjust=None, position='both', label=True, legend = True, fontsize=16, figsize=[7., 5.], tscale=1.,  verbose=True):
        '''
        the emPlot using the object oriented version of matplotlib - a figure and axis objects are returned
        to plot line intensities divided by gofnt
        adjust is to provide an adjustment to the position of the labels
        position : one of 'both', 'right', 'left', or 'none'
        this uses the modern object interface fig, ax = plt.subplots(figsize=figsize)
        fig, ax are in returned in the self.EmPlotOjb dict

        Keyword Arguments
        -----------------

        vs: `str`, either 'T', or 'D'
            whether to plot the emission measure vs temperature or density

        loc:  `str`
            matplotlib argument for plt.legend

        fs:  `int`
            the fontsize for the legend

        adjust:  `list`
            a list of multiplicative adjustments to the labels to the plot lines
            must be the same length as the number of lines

        position:  `str`
            where the labels to the lines should be placed, `both` for both ends, `left` for the left size only, 'right' for the right side only, or None for no labels

        label:  `bool`
            whether to apply

        legend:  `bool`
            whether to include a matplotlib legend

        fontsize:  `int`
            fontsize for the matplotlib xlabel and ylabel

        figsize:  two element `list` or `ndarray`
            sets the figure size when using matplotlib subplots to initiate the object style plotting

        tscale:  `float`
            scale the temperature by dividing by tscale

        verbose : `bool`
            if True, additional output is sent to the terminal

        '''
        match = self.Match
        temp = self.Temperature
        dens = self.EDensity
        nInt = len(match)
        if position == 'both':
            hpos = ['b']*nInt
        elif position == 'right':
            hpos = ['r']*nInt
        elif position == 'left':
            hpos = ['l']*nInt
        elif position is None:
            hpos = ['n']*nInt
        elif type(position) is list:
            if len(position) == nInt:
                hpos = position
                print(' position is the right length')
            else:
                print(' position is not the right length')
                return
        else:
            print(' position not understood: should be left, right, both, None')
            return

        if adjust is None:
            adjust = np.ones(nInt, np.float64)
    #        print(' nInt = %5i'%(nInt))
        if not 'intensity' in match[0]:
            print(' must run mgofnt or gofnt first')
            return
        elif vs == 'T':
            fig,  ax = plt.subplots(figsize=figsize)
            if temp[0] == temp[-1]:
                ntemp = 1
            else:
                ntemp = temp.size
#                ntemp = match[0]['intensity'].shape[0]
            if ntemp > 1:
                em = np.zeros((nInt, ntemp), 'float64')
                xvar = temp/tscale
            for idx in range(nInt):
                nonzed = match[idx]['intensitySum'] > 0.
                large = match[idx]['intensitySum'] >  match[idx]['intensitySum'].max()*0.01
                realgood = np.logical_and(nonzed, large)
                if realgood.sum() > 0 and match[idx]['obsIntensity'] > 0.:
                    wvlstr = ' %5.1f'%(match[idx]['obsWvl'])
                    em[idx][realgood] = match[idx]['obsIntensity']/match[idx]['intensitySum'][realgood]
                    if label:
                        ax.loglog(xvar[realgood], em[idx][realgood], lw=2, label = wvlstr)
                    else:
                        ax.loglog(xvar[realgood], em[idx][realgood], lw=2)
                    if ntemp > 1:
                        lblx1 = temp[realgood][0]
                        lblx2 = temp[realgood][-1]
                    else:
                        lblx1 = self.Temperature[realgood].min()
                        lbly1 = em[idx][0]
                    lbly1 = em[idx][realgood][0]
                    lbly2 = em[idx][realgood][-1]
                    if hpos[idx] == 'l':
                        ax.text(lblx1, lbly1, wvlstr.strip())
                    elif hpos[idx] == 'r':
                        ax.text(lblx2, lbly2, wvlstr.strip())
                    elif hpos[idx] == 'b':
                        ax.text(lblx2, lbly2, wvlstr.strip())
                        ax.text(lblx1, lbly1, wvlstr.strip())
                else:
                    print(' len of match[idx]wvl   %i'%(len(match[idx]['wvl'])))
                    if len(match[idx]['wvl']) != 0:
                        print('no values for idx = %5i wvl = %8.2f'%(idx, match[idx]['wvl']))
                        print(' nonzed = %i'%(nonzed.sum()))
                        print('intensity = %10.2e'%(match[idx]['intensity']))
            if legend:
                ax.legend(loc=loc, fontsize=fontsize)
            ax.set_ylabel('Emission Measure (cm$^{-5}$)', fontsize=fontsize)
            ax.set_xlabel('Temperature (K)',  fontsize=fontsize)
        elif vs != 'T':
            fig,  ax = plt.subplots(figsize=figsize)
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
                        ax.loglog(xvar[realgood], em[idx][realgood], lw=2, label=wvlstr)
                    else:
                        ax.loglog(xvar[realgood], em[idx][realgood], lw=2)
                    print(' %i  %5.1f'%(idx, match[idx]['obsWvl']))
                    lblx = dens
                    lbly = em[idx][realgood]*1.05
                    if hpos[idx] =='b':
                        ax.text(lblx[0], lbly[0]*adjust[idx], wvlstr.strip(), fontsize=14)
                        ax.text(lblx[-1], lbly[-1]*adjust[idx], wvlstr.strip(), horizontalalignment='right', fontsize=14)
                    elif hpos[idx] == 'r':
                        ax.text(lblx[-1], lbly[-1]*adjust[idx], wvlstr.strip(), horizontalalignment='right', fontsize=14)
                    elif hpos[idx] == 'l':
                        ax.text(lblx[0], lbly[0]*adjust[idx], wvlstr.strip(), fontsize=14)

                else:
                    print('no values for idx = %5i wvl = %8.2f'%(idx, match[idx]['wvl']))
                    print(' nonzed = %i'%(nonzed.sum()))
                    print('intensity = %10.2e'%(match[idx]['obsIntensity']))
            if legend:
                ax.legend(loc=loc, fontsize=fs)
            ax.set_ylabel('Emission Measure (cm$^{-5}$)', fontsize=fontsize)
            ax.set_xlabel('Electron Density (cm$^{-3}$)',  fontsize=fontsize)
            fig.tight_layout()
        self.EmPlotObj = {'fig':fig, 'ax':ax}
#        return fig,  ax
        #
        # --------------------------------------------------------------------------
        #
    def diff(self, sort=None,  verbose=False):
        '''
        calculates the weighted and straight differences between observed and predicted
        creates an  attribute self.Dict, a dict with the following keys:
        'wvl' = observed wavelength (A)
        'relDiff' = (I_obs - I_pred)/(I_obs)
        'ionS' the CHIANTI type name for an ion
        sort be either of none, 'wvl', or 'ion'


        Keyword Arguments
        -----------------

        sort:  `str` or None
            whether the output should be sorted by `wvl` or `ion` or not


        '''
        wghtFactor = self.WghtFactor
        nMatch = len(self.Match)
        if sort is None:
            sorter = range(nMatch)
        elif sort == 'ion':
            indexer = []
            for amatch in self.Match:
                ionStr = amatch['exptIon']
                ionDict = util.convertName(ionStr)
                Z = ionDict['Z']
                stage = ionDict['Ion']
                number = Z*1000 + stage
                indexer.append(number)
            sorter = np.argsort(indexer)
        elif sort == 'wvl':
            indexer = []
            for amatch in self.Match:
                wvlObs = amatch['obsWvl']
                indexer.append(wvlObs)
            sorter = np.argsort(indexer)

        wvlDiff = []
        diff = []
        wDiff = []
        intOverPred = []
        diffOverInt = []
        predOverInt = []

        noPredIdx = []
        noPredWvl = []
        noPredIon = []



        for iwvl in sorter:
            amatch = self.Match[iwvl]
            if amatch['predicted'] > 0. :
                wvlDiff.append(amatch['obsWvl'])
                diff.append(self.Intensity[iwvl]-amatch['predicted'])
                chi = np.abs(self.Intensity[iwvl]-amatch['predicted'])/(wghtFactor*self.Intensity[iwvl])
                wDiff.append(chi)
                intOverPred.append(self.Intensity[iwvl]/amatch['predicted'])
                predOverInt.append(amatch['predicted']/self.Intensity[iwvl])
                diffOverInt.append((self.Intensity[iwvl]-amatch['predicted'])/self.Intensity[iwvl])

            else:
                if verbose:
                    print(' no prediction wvl:  %8.2f ion:  %s'%(amatch['obsWvl'], amatch['ionS']))
                noPredIdx.append(iwvl)
                noPredWvl.append(self.Match[iwvl]['wvl'])
                noPredIon.append(self.Match[iwvl]['ions'])

        diffNp = np.asarray(diff, np.float64)
        intOverPredNp = np.asarray(intOverPred, np.float64)
        diffOverIntNp = np.asarray(diffOverInt, np.float64)
        predOverIntNp = np.asarray(predOverInt, np.float64)

        threeSig = 3.*diffOverIntNp.std()

        poor = np.abs(diffOverIntNp) > threeSig

        self.Diff = {'diff':diffNp, 'intOverPred':intOverPredNp, 'diffOverInt':diffOverIntNp,
            'wvl':wvlDiff, 'ionS':self.IonS, '3sig':threeSig, 'poor':poor, 'noPredIdx':noPredIdx,
            'noPredWvl':noPredWvl, 'noPredIon':noPredIon, 'predOverInt':predOverIntNp}
        #
        # -------------------------------------------------------------------------
        #
    def diffPlot(self, ratio = 'diffOverInt',  title=False, loc='upper right',  fontsize=16, figsize=[7., 5.]):
        """

        Parameters
        ----------

        ratio:  str
            determines which Y values to print
            can be 'diffOverInt' (default), or 'intOverPred', or 'predOverInt'

        title:  bool
            whether to plot the title or not

        fontsize:  int
            fontsize for matplotlib plots

        figsize:  2d list, ndarray
            the figure size for the plot


        Attributes
        ----------

        DiffPlot:  dict
            contains the fig, ax matplotlib objects created

        """
        goodRatio = ['diffOverInt',  'intOverPred', 'predOverInt']
        ylabels = ['(Obs - Pred)/Obs',  'Int / Pred', 'Pred / Int']
        if ratio not in goodRatio:
            print(' ratio = %s'%(ratio))
            print('should be set to %s  or  %s or %s'%(goodRatio[0],  goodRatio[1],  goodRatio[2]))
            return
        else:
            ratioIndex = goodRatio.index(ratio)
        if hasattr(self, 'Diff'):
            wvl = self.Diff['wvl']
#            diff = self.Diff['diffOverInt']
            diff = self.Diff[ratio]
            diffMean = diff.mean()
            diffStd = diff.std()
            fig,  ax = plt.subplots(figsize=figsize)
            ax.plot(wvl, diff,'o')
            ax.axhline(diffMean, color='k', lw=2, label='Mean')
            ax.axhline(diffMean + diffStd, color='r', lw=2, linestyle='--', label='1 std')
            ax.axhline(diffMean - diffStd, color='r', lw=2, linestyle='--')  #, label='1 std')
            ax.axhline(diffMean + 2.*diffStd, color='b', lw=2, linestyle='dotted', label='2 std')
            ax.axhline(diffMean - 2.*diffStd, color='b', lw=2, linestyle='dotted')  #, label='2 std')
            ax.axhline(diffMean + 3.*diffStd, color='g', lw=2, linestyle='dotted', label='3 std')
            ax.axhline(diffMean - 3.*diffStd, color='g', lw=2, linestyle='dotted')  #, label='3 std')
            ax.set_xlabel('Wavelength (\u212B)', fontsize=14)

            ax.set_ylabel(ylabels[ratioIndex], fontsize=14)
            if title:
                mytitle = 'diff Mean %10.3f  diff Std  %10.3f'%(diffMean, diffStd)
                ax.set_title(mytitle, fontsize=fontsize)
            ax.legend(loc=loc, fontsize=12) #  bbox_to_anchor=(0.99, 1.0),
            fig.tight_layout()
            self.DiffPlot = {'fig':fig, 'ax':ax}
        else:
            print(' the Diff attribute does not exist')
            print(' run method diff to create it')
        #
        # --------------------------------------------------------------------------
        #
    def diffPrint(self, filename='diffPrint.txt',  sort=None):
        '''
        calculates the weighted and straight differences between observed and predicted
        prints the values and saves to a file
        also created a attribute self.Dict, a dict with the following keys:
        'wvl' = observed wavelength (A)
        'relDiff' = (I_obs - I_pred)/(I_obs)
        'ionS' the CHIANTI type name for an ion


        Keyword Arguments
        -----------------

        filename:  `str`
            the filename where the text should be output

        sort:  `str`, can be `wvl`, `ion`, or None
            whether the output should be sorted by `wvl` or `ion` or not


        '''
        wghtFactor = self.WghtFactor
        cwd = os.getcwd()
        today = date.today()
        thisday =today.strftime('%Y_%B_%d')
        nMatch = len(self.Match)
        if sort is None:
            sorter = range(nMatch)
        elif sort == 'ion':
            indexer = []
            for amatch in self.Match:
                ionStr = amatch['exptIon']
                ionDict = util.convertName(ionStr)
                Z = ionDict['Z']
                stage = ionDict['Ion']
                number = Z*1000 + stage
                indexer.append(number)
            sorter = np.argsort(indexer)
        elif sort == 'wvl':
            indexer = []
            for amatch in self.Match:
                wvlObs = amatch['obsWvl']
                indexer.append(wvlObs)
            sorter = np.argsort(indexer)

        wvlDiff = []
        diff = []
        wDiff = []
        intOverPred = []
        diffOverInt = []
        predOverInt = []

        noPredIdx = []
        noPredWvl = []
        noPredIon = []

        dash = ' -------------------------------------------------'
        pformat1 = ' %5i %7s %10.3f %10.2e %10.2e %10.3f %10.3f %10.3f %10.3f'
        pformat1a = ' %5i %7s %10.3f %10.2e %10.2e *****NID******'
        sformat  = ' %5s %7s %10s %10s %10s %10s %10s %10s %10s'
        fsformat  = ' %5s %7s %10s %10s %10s %10s %10s %10s %10s\n'
        with open(filename, 'w') as outpt:
            if hasattr(self,  'SearchData'):
                idx = self.SearchData['best']['idx']
                dens = self.SearchData['best']['density']
                temp = self.SearchData['best']['temperature']
                emfit = self.SearchData['best']['emfit']
                em = self.SearchData['best']['em']
                print(' results from search')
                print(' %5s %10s %10s %10s %10s'%('index',  'density', 'temp', 'emfit','em'))
                pformat = '%5i %10.2e %10.2e %10.3f %10.2e'
                wformat = '%5i %10.2e %10.2e %10.3f %10.2e \n'
                outpt.write('results from search \n')
                outpt.write(' %5s %10s %10s %10s %10s \n'%('index',  'density', 'temp', 'emfit','em'))
                if type(idx) is int:
                    # this is for the case of a 1d search
                    idat = 0
                    index = idx
                    print(pformat%(index, dens, temp,  emfit,  em))
                    outpt.write(wformat%(index, dens, temp,  emfit,  em))
                else:
                    for idat, index in enumerate(idx):
                        print(pformat%(index, dens[idat], temp[idat],  emfit[idat],  em[idat]))
                        outpt.write(wformat%(index, dens[idat], temp[idat],  emfit[idat],  em[idat]))
#                print(' %5i %10.2e %10.2e %10.3f %10.2e'%(idx, dens, temp, emfit, em))
#                outpt.write(' %5i %10.2e %10.2e %10.3f %10.2e \n'%(idx, dens, temp, emfit, em))
            else:
                emIndices = self.EmIndices
                for emidx in emIndices:
                    print(' %5i %12.2e %12.2e %12.2e'%(emidx, self.EDensity[emidx],  self.Temperature[emidx], self.Em[emidx]))
                    outpt.write(' %5i %12.2e %12.2e %12.2e \n'%(emidx, self.EDensity[emidx],  self.Temperature[emidx], self.Em[emidx]))
            print(dash)
            outpt.write(dash  + '\n')
            print(' cwd:  %s'%(cwd))
            outpt.write(' cwd:  %s \n'%(cwd))
            print(' today is %s'%(thisday))
            outpt.write(' today is %s \n'%(thisday))
            print(' WghtFactor = %10.3f'%(wghtFactor))
            outpt.write(' WghtFactor = %10.3f \n'%(wghtFactor))
            emIndices = self.EmIndices
            print(' %5s %12s %12s %12s'%('index',  'density',  'temperature',  'Em'))
            outpt.write(' %5s %12s %12s %12s \n'%('index',  'density',  'temperature',  'Em'))
            for emidx in emIndices:
                print(' %5i %12.2e %12.2e %12.2e %8.3f'%(emidx, self.EDensity[emidx],  self.Temperature[emidx], self.Em[emidx],  self.EmLog[emidx]))
                outpt.write(' %5i %12.2e %12.2e %12.2e %8.3f \n'%(emidx, self.EDensity[emidx],  self.Temperature[emidx], self.Em[emidx],  self.EmLog[emidx]))
            print(dash)
            outpt.write(dash+'\n')
            print(' chi = abs(int - pred)/(wght*int))  strDiff = (int - pred)/pred')
            outpt.write(' chi = abs(int-pred)/(wght*int))  strDiff = (int - pred)/int \n')
            print(sformat%('', '', 'A', '', '', '', '', 'abs', 'abs' ))
            print(sformat%('iwvl', 'ionS', 'wvl', 'intensity', 'predicted', 'int/pred', 'chi', 'relDev', 'dif/int'))
            outpt.write(fsformat%('', '', 'A', '', '', '', '', ' ', ' ' ))
            outpt.write(fsformat%('iwvl', 'ionS', 'wvl', 'intensity', 'predicted', 'int/pred', 'chi', 'relDev',  'dif/int'))
            for iwvl in sorter:
                amatch = self.Match[iwvl]
                if amatch['predicted'] > 0. :
                    wvlDiff.append(amatch['obsWvl'])
                    diff.append(self.Intensity[iwvl]-amatch['predicted'])
                    chi = np.abs(self.Intensity[iwvl]-amatch['predicted'])/(wghtFactor*self.Intensity[iwvl])
                    wDiff.append(chi)
                    intOverPred.append(self.Intensity[iwvl]/amatch['predicted'])
                    predOverInt.append(amatch['predicted']/self.Intensity[iwvl])
                    diffOverInt.append((self.Intensity[iwvl]-amatch['predicted'])/self.Intensity[iwvl])
                    pstring = pformat1%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'], self.Intensity[iwvl]/amatch['predicted'], chi, intOverPred[-1],  diffOverInt[-1] )
                else:
                    pstring = pformat1a%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'])
                    noPredIdx.append(iwvl)
                    noPredWvl.append(self.Match[iwvl]['wvl'])
                    noPredIon.append(self.Match[iwvl]['ions'])

                print(pstring)
                outpt.write(pstring+'\n')

            print(' WghtFactor = %10.3f'%(wghtFactor))
            outpt.write(' WghtFactor = %10.3f \n'%(wghtFactor))

            diffNp = np.asarray(diff, np.float64)
            intOverPredNp = np.asarray(intOverPred, np.float64)
            predOverIntNp = np.asarray(predOverInt, np.float64)
            diffOverIntNp = np.asarray(diffOverInt, np.float64)
            wghtDiffOverIntNp = diffOverIntNp/wghtFactor
            chisq2 = (wghtDiffOverIntNp**2).sum()

            print(' mean of relative Deviation = %10.3f 3*relDev = %10.3f stdDev = %10.3f\n'%(intOverPredNp.mean(), 3.*intOverPredNp.mean(), intOverPredNp.std()))
            outpt.write(' mean of rel Dev = %10.3f 3*rel Dev  = %10.3f stdDiff:  %10.3f \n'%(intOverPredNp.mean(), 3.*intOverPredNp.mean(), intOverPredNp.std()))
            print(' mean of intOverPred = abs(int/pred -1.) %10.3f'%(intOverPredNp.mean()))
            outpt.write(' mean of intOverPred = abs(int/pred -1.) %10.3f \n'%(intOverPredNp.mean()))

            onestd = diffOverIntNp.std()
            threeSig = 3.*diffOverIntNp.std()
            print(' std = %10.3f 3*std = %10.3f'%(onestd,  threeSig))
            outpt.write(' std = %10.3f   3*std = %10.3f \n'%(onestd,  threeSig))

            print(' sum of diffOverInt^2 =  %10.3f weighted %10.3f '%((diffOverIntNp**2).sum(),  (diffOverIntNp**2).sum()/wghtFactor**2))
            outpt.write(' sum of diffOverInt^2 =  %10.3f  weighted %10.3f\n'%((diffOverIntNp**2).sum(), (diffOverIntNp**2).sum()/wghtFactor**2  ))

            print(' chisq2 = %10.3f'%(chisq2))
            outpt.write(' chisq2 = %10.3f \n'%(chisq2))

            print(' Nobs = %i Nions = %i'%(self.Nobs,  self.Nions))
            outpt.write(' Nobs = %i  Nions = %i \n'%(self.Nobs,  self.Nions))
            print(' Nparams = %i'%(self.Nparams))
            outpt.write(' Nfree = %i \n'%(self.Nparams))
            chisq,  msk = self.getChisq()
            normChisq = self.getNormalizedChisq()
            print('           Chisq = %10.3f'%(chisq))
            print('Normalized Chisq = %10.3f chisq/(nobs)'%(chisq/float(nMatch)))
            print('Reduced Chisq    = %10.3f chisq/(nobs - nparams)'%(chisq/(nMatch - self.Nparams)))
            outpt.write('           Chisq = %10.3f \n'%(chisq))
            outpt.write('Normalized Chisq = %10.3f chisq/(nobs) \n'%(normChisq))
            outpt.write('Reduced Chisq    = %10.3f chisq/(nobs - nparams)\n'%(chisq/(nMatch - self.Nparams)))
            poor = np.abs(diffOverIntNp) > threeSig
#            idx = np.arange(nMatch)
#            pdx = idx[poor]
            if poor.size != len(sorter):
                print(' poor.size %i  len(sorter) %i'%(poor.size,  len(sorter)))
            npsorter = np.asarray(sorter)
            pdx = npsorter[poor]
            for i in pdx:
                print('%5i %s %10.3f %10.3f %10.3f'%(i, self.IonS[i],  self.Wvl[i], np.abs(intOverPredNp)[i],  diffOverIntNp[i]))
                outpt.write('%5i %s %10.3f %10.3f %10.3f\n'%(i, self.IonS[i],  self.Wvl[i], np.abs(intOverPredNp)[i], diffOverIntNp[i]))
        self.Diff = {'diff':diffNp, 'intOverPred':intOverPredNp, 'diffOverInt':diffOverIntNp,
            'wvl':wvlDiff, 'ionS':self.IonS, '3sig':threeSig, 'poor':poor, 'noPredIdx':noPredIdx,
            'noPredWvl':noPredWvl, 'noPredIon':noPredIon, 'predOverInt':predOverIntNp}
        #
        # --------------------------------------------------------------------------
        #
    def predict(self):
        '''
        to predict the intensities of the observed lines from an emission measure
        the emission measure is already specified as self.Em which is an ndarray
        the temperatures are set by emSetIndices
        '''
        #
        for iwvl, amatch in enumerate(self.Match):
            pred = np.sum(amatch['intensitySum']*self.Em)
            try:
                self.Match[iwvl]['predicted'] = pred
            except:
                self.Match[iwvl]['predicted'] = 0.
                print(' error in predict, iwvl = %5i'%(iwvl))
        #
        # --------------------------------------------------------------------------
        #
    def predictPrint(self, minContribution=0.1, filename='predictPrint.txt', sort=None, verbose=0):
        '''
        to predict the intensities of the observed lines from an emission measure
        the emission measure is already specified as self.Em which is an np array
        sort can be 'wvl' or 'ion', otherwise, there is no sorting done


        Keyword Arguments
        -----------------

        minContribution:  `float`
            the minimum contribution a blend must supply to be included in the text output

        filename:  `str`
            the filename where the text should be output

        sort:  `str`, can be `wvl`, `ion`, or None
            whether the output should be sorted by `wvl` or `ion` or not

        verbose : `bool`
            if True, additional output is sent to the terminal

        '''
        cwd = os.getcwd()
        wghtFactor = self.WghtFactor
        nMatch = len(self.Match)
        if verbose:
            print('nMatch:  %i'%(nMatch))
        if sort is None:
            sorter = range(nMatch)
        elif sort == 'ion':
            indexer = []
            for amatch in self.Match:
                ionStr = amatch['exptIon']
                ionDict = util.convertName(ionStr)
                Z = ionDict['Z']
                stage = ionDict['Ion']
                number = Z*1000 + stage
                indexer.append(number)
            sorter = np.argsort(indexer)
        elif sort == 'wvl':
            indexer = []
            for amatch in self.Match:
                wvlObs = amatch['obsWvl']
                indexer.append(wvlObs)
            sorter = np.argsort(indexer)

        dash = ' -------------------------------------------------'
        # for cases where no predicted intensity
        pformatNone = ' %5i %7s %10.3f'
        pformat1 = ' %5i %7s %10.3f %10.2e %10.2e %10.3f %10.3f'
        pformat1s = ' %5s %7s %10s %10s %10s %10s %10s'
        pformat2 = '         %s'
        pformat3 = '        %10.3f %4i %4i %20s - %20s %5i %5i %7.3f'
        pformat3s = '        %10s %4s %4s %20s - %20s %5s %5s %s'
        if filename:
            with open(filename, 'w') as outpt:
                print(' cwd:  %s'%(cwd))
                outpt.write(' cwd:  %s \n'%(cwd))
                try:
                    emIndices = self.EmIndices
                    for emidx in emIndices:
                        print(' %5i %12.2e %12.2e %12.2e'%(emidx, self.EDensity[emidx],  self.Temperature[emidx], self.Em[emidx]))
                        outpt.write(' %5i %12.2e %12.2e %12.2e\n'%(emidx, self.EDensity[emidx],  self.Temperature[emidx], self.Em[emidx]))
                except:
                    print(' attribute EmIndices has not been created')
                    outpt.write(' attribute EmIndices has not been created \n')
                    return
                if hasattr(self, 'MatchName'):
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
                    amatch = self.Match[iwvl]
                    if amatch['predicted'] > 0. :
                        chi = (np.abs(self.Intensity[iwvl]-amatch['predicted'])/(wghtFactor*self.Intensity[iwvl]))
                        pstring = pformat1%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'], self.Intensity[iwvl]/amatch['predicted'], chi  )
                    else:

                        pstring = pformat1%(iwvl, self.IonS[iwvl],  self.Wvl[iwvl], self.Intensity[iwvl], amatch['predicted'], -1.,  -1.)
                    print(pstring)
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
                                    print(pformat2%(anion))
                                    outpt.write(pformat2%(anion)+'\n')
                                    ionPrint = 1
                                print(pformat3%(awvl, amatch['lvl1'][jon][iline], amatch['lvl2'][jon][iline], amatch['pretty1'][jon][iline], amatch['pretty2'][jon][iline].ljust(20), amatch['lineIdx'][jon][iline], amatch['predictedLine'][jon][iline], contrib))
                                outpt.write(pformat3%(awvl, amatch['lvl1'][jon][iline], amatch['lvl2'][jon][iline], amatch['pretty1'][jon][iline], amatch['pretty2'][jon][iline].ljust(20), amatch['lineIdx'][jon][iline], amatch['predictedLine'][jon][iline], contrib)+'\n')
                                print(dash)
                                outpt.write(dash +'\n')
                    print(dash)
                outpt.write(dash +'\n')
        #
        # --------------------------------------------------------------------------
        #
    def predictPrint1d(self, minContribution=0.1, filename='predictPrint1d.txt', verbose=False):
        '''
        to predict the intensities of the observed lines from an emission measure
        the emission measure is already specified as self.Em which is an np array

        to be used after a 1d search over density


        Keyword Arguments
        -----------------

        minContribution:  `float`
            the minimum contribution a blend must supply to be included in the text output

        filename:  `str`
            the filename where the text should be output

        verbose : `bool`
            if True, additional output is sent to the terminal


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
        for iwvl,  amatch in enumerate(self.Match):
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
        if filename:
            with open(filename, 'w') as outpt:
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
                for iwvl,  amatch in enumerate(self.Match):
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
        return the weighted chi-squared
        '''
        weightedDiff,  msk = self.getWeightedDiff()
        chisq = np.sum(weightedDiff**2)
        return chisq, msk

    def getNormalizedChisq(self):
        '''
        return normalized chisq:  chi-squared divided by the number of observed lines
        '''
        chisq,  mask = self.getChisq()
        normalChisq = chisq/float(self.Nobs)
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
        nwvl = len(self.Match)
        weightedDiff = np.zeros(nwvl, 'float64')
        for iwvl, amatch in enumerate(self.Match):
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
        greater than 0


        Keyword Arguments
        -----------------

        verbose : `bool`
            if True, additional output is sent to the terminal

        '''
        nT = len(self.Temperature)
        nlines = len(self.Match)
        print(' n lines = %5i '%(nlines))
        minDx = 0
        maxDx = nT
        for iline, match in enumerate(self.Match):
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
        called by fit1t

        Parameters
        ----------

        em:  number
            the log10 value of the emission measure

        Returns
        -------

        weighted chisquared:  float


        Todo
        ----

        see if this can be replaced by fitFuncNt


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
        the fitting function for the multiple  temperature model to be called by leastsq
        called by fitNt

        Parameters
        ----------

            value:  `list`
                the initial values for the em fit
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
        used in search1dspace

        Parameters
        ----------

        initialValue: `float`
            the initial value to start the leastsq process


        Todo
        ----

        see if this can be replaced by fitFunct1t or fitNt
        '''
        out = optimize.leastsq(self.fitFunc1t, np.asarray(initialValue, 'float64'), full_output=1, maxfev=maxfev)
        self.Leastsq = {'em':out[0], 'cov':out[1], 'info':out[2], 'message':out[3], 'ier':out[4]}
        #
        # -------------------------------------------------------------
        #
    def fitNt(self, initialValue, maxfev=0):
        '''
        calls leastsq to fit the multi temperature models
        called by search2tSpace, search3tSpace etc

        Parameters
        ----------

        initialValue:  `list`
            the initial trial value for the emission measure (log1))

        Keyword Arguments
        -----------------

            maxfev:  `int`
                not sure it is needed
        '''
        out = optimize.leastsq(self.fitFuncNt, np.asarray(initialValue, 'float64'), full_output=1, maxfev=maxfev)
        self.Leastsq = {'em':out[0], 'cov':out[1], 'info':out[2], 'message':out[3], 'ier':out[4]}
        #
        # -----------------------------------------------------------
        #
    def search1dSpace(self, initialEm, indxlimits=None, verbose=False, log=False, maxfev=0):
        '''
        to conduct a brute force search over electron density for an isothermal-space and find the
        best fit to the em and density
        indxlimits give the range of indices to fit over
        can use self.MinIndex and self.MaxIndex+1
        initialEm = log value of the emission measure to begin the searching

        Parameters
        ----------

        initialEm:  `float`
            the initial trial value for the log10 emission measure

        Keyword Arguments
        -----------------

        indxlimits:  `list`, None
            the range of indices of the density array to search
            if None is specified, the whole range is searched

        verbose:  `bool`
            if True, additional output is sent to the terminal

        log:  `bool`
            if True, a log file is created - 'search1d.raw'

        '''
        t1 = datetime.now()
        if self.Nobs <= 2.:
            print(' the number of observables %i is less than the number of parameters'%(self.Nobs))
        if indxlimits is None:
            if not hasattr(self, 'MinIndex'):
                self.findMinMaxIndices()
            indxlimits = [self.MinIndex, self.MaxIndex]
        self.Nparams = 2. + 1.

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
        reducedChisq = minChisq/float(self.Nobs - self.Nparams)
        #
        # key 'emfit' is the log value, 'em' is actual value
        self.SearchData['best'] = {'em':em, 'emfit':emfit[gdx[0]], 'chisq':chisq[gdx[0]], 'reducedChisq':reducedChisq, 'idx':idx[gdx[0]], 'density':self.EDensity[idx[gdx[0]]], 'temperature':self.Temperature[idx[gdx[0]]]}
        t2 = datetime.now()
        dt = t2 - t1
        print(' elapsed seconds = %12.3f'%(dt.seconds))
        #
        #-----------------------------------------------------
        #
    def search1tEmSpace(self, verbose=False):
        ''' to find the value of chisq as a function of Em with T = best-fit

        Keyword Arguments
        -----------------

        verbose:  `bool`
            if True, additional output is sent to the terminal

        '''
        self.Nparams = 2. + 1.
        if not hasattr(self, 'SearchData'):
            print(' dem has not been searched yet')
            return
        if self.Nobs <= 2.:
            print(' the number of observables %10.1f is less than the number of parameters %10.1f'%(self.Nobs, 2.*2.))
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
        self.SearchData['EmError'] = {'em':emSpace, 'emfit':emfitSpace, 'chisq':emChisq, 'idx':self.SearchData['best']['idx'], 'temperature':self.Temperature[self.SearchData['best']['idx']], 'nparams':self.Nparams}
        #
        #-----------------------------------------------------
        #
    def search2tSpace(self, initial, indxlimits=None, verbose=0, log=0, maxfev=0):
        '''
        to conduct a brute force search of 2 temperature space and find the
        best fit
        indxlimits give the range of indices to fit over


        Parameters
        ----------

        initial:  `list`
            the initial trial values (2) for the log10 emission measure

        Keyword Arguments
        -----------------

        indxlimits:  `list`, None
            the range of indices of the density array to search
            if None is specified, the whole range is searched

        verbose:  `bool`
            if True, additional output is sent to the terminal

        log:  `bool`
            if True, a log file is created - 'search1d.raw'

        '''
        self.Nparams = 4. + 1.
        t1 = datetime.now()
        if self.Nobs <= 2.*2.:
            print(' the number of observables %10.2e is less than the number of parameters %10.1f'%(self.Nobs, 2.*2.))
            return
        if indxlimits == None:
            if not hasattr(self, 'MinIndex'):
                self.findMinMaxIndices()
            indxlimits = [self.MinIndex, self.MaxIndex]
        print('indxlimits = %5i %5i '%(indxlimits[0], indxlimits[1]))
        nTemp = len(self.EmIndices)
        if log:
            logname = 'search2t.raw'
            logfile = open(logname, 'w')
            logfile.write('%s \n'%(self.SpecData['filename']))
#        pformat = '%5i %4i %4i %4i %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e'
        logformat = '%5i %4i %4i %4i %4i %10.2e %10.2e %10.2e %10.2e %10.2e \n'

        temp1Searched = self.Temperature[indxlimits[0]:indxlimits[1]]
        temp2Searched = self.Temperature[indxlimits[0]+1:indxlimits[1]+1]

        print('temp1Searched.size, temp2Searched.size %5i %5i'%(temp1Searched.size, temp2Searched.size))

        emfit = []
        idx12 = []
        idx = []
        searchDx12 = []
        chisq = []
        info = []
        maskedValues = []
#        msk = np.ma.make_mask(np.ones((temp1Searched.size, temp2Searched.size)))
#        msk2 = np.ma.make_mask(np.ones((temp1Searched.size, temp2Searched.size, 2)))

        chisq2d = np.ma.array(np.zeros((temp1Searched.size, temp2Searched.size), 'float64'), mask=True)
        temp2d = np.zeros((nTemp, temp1Searched.size, temp2Searched.size), 'float64')
        em2d = np.ma.array(np.zeros((nTemp, temp1Searched.size, temp2Searched.size), 'float64'), mask=True)
#        emNd = np.ma.array(np.zeros((nTemp, temp1Searched.size, temp2Searched.size, temp3Searched.size), 'float64'), mask=True)

        counter = 0
        for idx1 in range(indxlimits[0], indxlimits[1], 1):
            for idx2 in range(idx1+1, indxlimits[1]+1, 1):
                idx12.append([idx1, idx2])
                idx.append([idx1, idx2])
                searchDx1 = idx1 - indxlimits[0]
                searchDx2 = idx2 - indxlimits[0] - 1
                self.emSetIndices([idx1, idx2])
                self.fitNt(initial, maxfev=maxfev)
                self.emSetIndices([idx1, idx2])
                self.fitNt(initial, maxfev=maxfev)
                # if there is no solution for one of the em values, it is set to a negative no.
#                gtz = self.Leastsq['em'] > 0.
                chisq1,  msk1 = self.getChisq()
                chisq.append(chisq1)
                temp2d[0, searchDx1, searchDx2] = temp1Searched[searchDx1]
                temp2d[1, searchDx1, searchDx2] = temp2Searched[searchDx2]
                searchDx12.append([searchDx1, searchDx2])
                if not msk1:
                    emfit.append(self.Leastsq['em'])
                    chisq2d.data[searchDx1, searchDx2] = chisq1
                    chisq2d.mask[searchDx1, searchDx2] = False
                    for it in range(nTemp):
                        em2d.data[it, searchDx1, searchDx2] = self.Leastsq['em'][it]
                        em2d.mask[it, searchDx1, searchDx2] = False
                    info.append(self.Leastsq['info'])
                    if verbose:
                        print(' %5i %5i %5i %5i %5i %12.2e %10.3f %12.2e %10.3f %12.2e %10i'%(counter, idx1, idx2, searchDx1, searchDx2, self.Temperature[idx1], self.Leastsq['em'][0], self.Temperature[idx2], self.Leastsq['em'][1], chisq[-1], self.Leastsq['info']['nfev']))
                    elif log:
                        logfile.write(logformat%(counter, idx1, idx2, searchDx1, searchDx2, self.Temperature[idx1], self.Leastsq['em'][0], self.Temperature[idx2], self.Leastsq['em'][1], chisq[-1]))
                else:
                    maskedValues.append([idx1, idx2])
                    if verbose:
                        print(' %5i %5i %5i bad em values %10.2f %10.2f'%(counter, idx1, idx2, self.Leastsq['em'][0], self.Leastsq['em'][1]))
                    elif log:
                        logfile.write('msk at %5i %4i %4i %10.2e \n'%(counter, idx1, idx2, chisq[-1]))
                    counter += 1
                    # these values already are masked
                    emfit.append(self.Leastsq['em'])
                    chisq2d.data[searchDx1, searchDx2] = chisq1
                    chisq2d.mask[searchDx1, searchDx2] = True
                    for it in range(nTemp):
                        em2d.data[it, searchDx1, searchDx2] = self.Leastsq['em'][it]
                        em2d.mask[it, searchDx1, searchDx2] = True
                counter += 1
        logfile.close()
        #
        self.SearchData = {'temperature':self.Temperature,'temp1Searched':temp1Searched,
            'temp2Searched':temp2Searched, 'emfit':emfit, 'em2d':em2d, 'idx12':idx12, 'idx':idx,
            'searchDx12':searchDx12, 'chisq':chisq, 'minchisq':min(chisq), 'chisq2d':chisq2d,
            'maskedValues':maskedValues, 'temp2d':temp2d,
            'message':'this contains all the data for the 2tExpSpace', 'nparams':self.Nparams}
        #
        print('min chisq = %10.3e at %5i '%(min(chisq), np.argmin(chisq)))
        print(' # of masked values %i'%(len(maskedValues)))
        #
        #  set things to best fit
        #
        gdx=[i for i,ch in enumerate(chisq) if ch == min(chisq)]
        if verbose:
#        print(' gdx = ', gdx)
            print(' len, min of chisq %5i %10.3e'%(len(chisq), min(chisq)))
            if len(gdx) == 1:
                print('single value gdx = %i '%(gdx[0]))  #, idx12[gdx[0]][0], idx12[gdx[0]][1])
#                temperature1 = eis.SearchData['temp1Searched']
            else:
                pstring = 'mult values of gdx = '
                for one in gdx:
                    pstring += '  %10i'%(one)
#        print(' emfit', emfit[gdx[0]])
#        print(' chisq', chisq[gdx[0]])
#        print(' idx ', idx12[gdx[0]])
#        print(' temperature', self.Temperature[idx12[gdx[0]]])
#        print(' density %10.2e'%(self.Density))
#        self.em2tSetIndices(idx12[gdx[0]])
#        self.em2texp(emfit[gdx[0]])
        print('indices of minimum after search =  %i %i'%(idx12[gdx[0]][0], idx12[gdx[0]][1] ))
        # gives 90% confidence for 4 parameters
#        chisq1sig = min(chisq)*np.ones_like(chisq2d, 'float64') + 7.779
#        good1sig = chisq2d < min(chisq) + 7.779
#        temp1sig1 = temp2d[good1sig, 0]
#        temp2sig1 = temp2d[good1sig, 1]
#        temp1sig = [temperatureSearched[good1sig].min(), temperatureSearched[good1sig].max()]
        self.emSetIndices(idx12[gdx[0]])
        self.emSet(emfit[gdx[0]])
        self.predict()
        em = 10.**emfit[gdx[0]]
        #

        # next line is needed for predict and predict print
        self.emSetIndices(idx12[gdx[0]])
        temperature1 = self.Temperature[idx12[gdx[0]][0]]
        temperature2 = self.Temperature[idx12[gdx[0]][1]]

        emfit1 = emfit[gdx[0]]
        adx = idx12[gdx[0]]
        minChisq = chisq[gdx[0]]
        reducedChisq = minChisq/float(self.Nobs - 2.*nTemp)


        print(' idx1 %5i idx2 %5i  '%(adx[0], adx[1]))
        print(' Temp =  %10.2e %10.2e '%(self.Temperature[adx[0]], self.Temperature[adx[1]]))
        print(' EM =    %10.2e %10.2e '%(em[0], em[1]))
        print(' EmFit = %10.3f %10.3f  '%(emfit1[0], emfit1[1]))

        t2 = datetime.now()
        dt = t2 - t1
        print(' elapsed seconds = %12.3f'%(dt.seconds))

        self.SearchData['best'] = {'em':em, 'emfit':emfit[gdx[0]], 'chisq':chisq[gdx[0]],
            'idx12':idx12[gdx[0]], 'idx':idx[gdx[0]],  'reducedChisq':reducedChisq,
            'temperature':[temperature1, temperature2], 'temperature1':temperature1,
            'temperature2':temperature2, 'density':self.Density[idx12[gdx[0]]],
            'gdx':gdx, 'searchDx12':searchDx12[gdx[0]]}
        #
        # -----------------------------------------------------------
        #
    def search3tSpace(self, initial, indxlimits=None, verbose=0, log=0):
        '''
        to conduct a brute force search of 3 temperature space and find the
        best fit

        Parameters
        ----------

        initial:  `list`
            the initial trial values (3) for the log10 emission measure


        Keyword Arguments
        -----------------

        indxlimits:  `list`, None
            the range of indices of the density array to search
            if None is specified, the whole range is searched

        verbose:  `bool`
            if True, additional output is sent to the terminal

        log:  `bool`
            if True, a log file is created - 'search1d.raw'

        '''
        self.Nparams = 6. + 1.
        t1=datetime.now()

        if indxlimits == None:
            if not hasattr(self, 'MinIndex'):
                self.findMinMaxIndices()
            indxlimits = [self.MinIndex, self.MaxIndex]
        print('indxlimits = %5i %5i '%(indxlimits[0], indxlimits[1]))
        nTemp = len(self.EmIndices)
        emfit = []
        idx123 = []
        idx = []
        searchDx123 = []
        chisq = []
        maskedValues = []
        temp1Searched = self.Temperature[indxlimits[0]:indxlimits[1]-1]
        temp2Searched = self.Temperature[indxlimits[0]+1:indxlimits[1]]
        temp3Searched = self.Temperature[indxlimits[0]+2:indxlimits[1]+1]

#        msk  = np.ma.make_mask(np.ones((temp1Searched.size, temp2Searched.size, temp3Searched.size)))
#        msk2 = np.ma.make_mask(np.ones((temp1Searched.size, temp2Searched.size, temp3Searched.size, 3)))
        chisqNd = np.ma.array(np.zeros((temp1Searched.size, temp2Searched.size, temp3Searched.size), 'float64'), mask=True)
        tempNd = np.array(np.zeros((nTemp, temp1Searched.size, temp2Searched.size, temp3Searched.size), 'float64'))
        emNd = np.ma.array(np.zeros((nTemp, temp1Searched.size, temp2Searched.size, temp3Searched.size), 'float64'), mask=True)

        if log:
            logname = 'search3t.raw'
            logfile = open(logname, 'w')
            logfile.write('%s \n'%(self.SpecData['filename']))
        # need to make this look like search2t
        pformat = '%5i %4i %4i %4i %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e'
        logformat = '%5i %4i %4i %4i %4i %4i %4i %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e \n'
        counter = 0
        for idx1 in range(indxlimits[0], indxlimits[1] - 1, 1):
            for idx2 in range(idx1+1, indxlimits[1], 1):
                for idx3 in range(idx2+1, indxlimits[1] + 1, 1):
                    self.emSetIndices([idx1, idx2, idx3])
                    self.fitNt(initial)
        # need to make this look like search 2t
                    chisq1,  msk1 = self.getChisq()
                    searchDx1 = idx1 - indxlimits[0]
                    searchDx2 = idx2 - indxlimits[0] - 1
                    searchDx3 = idx3 - indxlimits[0] - 2
                    searchDx123.append([searchDx1, searchDx2, searchDx3])
                    tempNd[0, searchDx1, searchDx2, searchDx3] = temp1Searched[searchDx1]
                    tempNd[1, searchDx1, searchDx2, searchDx3] = temp2Searched[searchDx2]
                    tempNd[2, searchDx1, searchDx2, searchDx3] = temp3Searched[searchDx3]
                    if not msk1:
                        emfit.append(self.Leastsq['em'])
                        chisq.append(chisq1)
                        idx123.append([idx1, idx2, idx3])
                        idx.append([idx1, idx2, idx3])
                        chisqNd.data[searchDx1, searchDx2, searchDx3] = chisq1
                        chisqNd.mask[searchDx1, searchDx2, searchDx3] = False
                        for it in range(nTemp):
                            emNd.data[it, searchDx1, searchDx2, searchDx3] = self.Leastsq['em'][it]
                            emNd.data[it, searchDx1, searchDx2, searchDx3] = self.Leastsq['em'][it]
                            emNd.data[it, searchDx1, searchDx2, searchDx3] = self.Leastsq['em'][it]
                            emNd.mask[it, searchDx1, searchDx2, searchDx3] = False

                        if verbose:
                            print(pformat%(counter, idx1, idx2, idx3,  self.Temperature[idx1], self.Leastsq['em'][0], self.Temperature[idx2], self.Leastsq['em'][1], self.Temperature[idx3], self.Leastsq['em'][2], chisq[-1]))
                        elif log:
                            logfile.write(logformat%(counter, idx1, idx2, idx3,  searchDx1, searchDx2, searchDx3, self.Temperature[idx1], self.Leastsq['em'][0], self.Temperature[idx2], self.Leastsq['em'][1], self.Temperature[idx3], self.Leastsq['em'][2], chisq[-1]))

                    else:
                        emfit.append(self.Leastsq['em'])
                        idx123.append([idx1, idx2, idx3])
                        idx.append([idx1, idx2, idx3])
                        chisq.append(chisq1)
                        maskedValues.append([idx1, idx2, idx3])
                        chisqNd.data[searchDx1, searchDx2, searchDx3] = chisq1
#                        chisq3d.mask[searchDx1, searchDx2, searchDx3] = False
                        for it in range(nTemp):
                            emNd.data[it, searchDx1, searchDx2, searchDx3] = self.Leastsq['em'][it]
                            emNd.data[it, searchDx1, searchDx2, searchDx3] = self.Leastsq['em'][it]
                            emNd.data[it, searchDx1, searchDx2, searchDx3] = self.Leastsq['em'][it]
                            emNd.mask[it, searchDx1, searchDx2, searchDx3] = True

                        if verbose:
                            print('msk at %5i %4i %4i %4i %10.2e'%(counter, idx1, idx2, idx3, chisq[-1]))
                        elif log:
                            logfile.write('msk at %5i %4i %4i %4i %10.2e \n'%(counter, idx1, idx2, idx3, chisq[-1]))
                    counter += 1
        if log:
            logfile.close()
        #
        self.SearchData = {'temperature':self.Temperature,'temp1Searched':temp1Searched,
            'temp2Searched':temp2Searched, 'temp3Searched':temp3Searched, 'emfit':emfit, 'idx':idx123,
            'searchDx123':searchDx123, 'chisq':chisq, 'chisqNd':chisqNd, 'minchiseq':min(chisq),
            'tempNd':tempNd, 'emNd':emNd, 'message':'this contains all the data for the 3tExpSpace',
            'nparams':self.Nparams}
        #
        #
        #  set things to best fit
        #
        print('min chisq = %10.3e at %i'%(min(chisq), np.argmin(chisq)))

        gdx=[i for i,ch in enumerate(chisq) if ch == min(chisq)]
        if verbose:
#        print(' gdx = ', gdx)
            print(' len, min of chisq %5i %10.3e'%(len(chisq), min(chisq)))
            if len(gdx) == 1:
                print('single value gdx = %i '%(gdx[0]))  #, idx12[gdx[0]][0], idx12[gdx[0]][1])
#                temperature1 = eis.SearchData['temp1Searched']
            else:
                pstring = 'mult values of gdx = '
                for one in gdx:
                    pstring += '  %10i'%(one)
        # 90% confidence for 6 parameters
#        chisq1sig = min(chisq)*np.ones_like(chisqNd, 'float64') + 10.645
#        good1sig = chisqNd < min(chisq) + 10.645

        self.emSetIndices(idx123[gdx[0]])
        self.emSet(emfit[gdx[0]])
        self.predict()
        em = 10.**emfit[gdx[0]]
        emfit1 = emfit[gdx[0]]
        adx = idx123[gdx[0]]
        minChisq = chisq[gdx[0]]
        reducedChisq = minChisq/float(self.Nobs - 2.*nTemp)

        print(' idx1 %5i idx2 %5i idx3 %5i '%(adx[0], adx[1], adx[2]))
        print(' Temp =  %10.2e %10.2e %10.2e '%(self.Temperature[adx[0]], self.Temperature[adx[1]], self.Temperature[adx[2]]))
        print(' EM =    %10.2e %10.2e %10.2e '%(em[0], em[1], em[2]))
        print(' EmFit = %10.3f %10.3f %10.3f '%(emfit1[0], emfit1[1], emfit1[2]))
        #
        t2 = datetime.now()
        dt = t2 - t1
        print(' elapsed seconds = %12.3f'%(dt.seconds))

        self.SearchData['best'] = {'em':em,'emfit':emfit[gdx[0]], 'chisq':chisq[gdx[0]],
        'reducedChisq':reducedChisq, 'idx':idx123[gdx[0]], 'gdx':gdx[0],
        'temperature':self.Temperature[idx123[gdx[0]]], 'density':self.Density}
        #
        # -----------------------------------------------------------
        #
    def search4tSpace(self, initial, indxlimits=None, verbose=0, log=0):
        '''
        to conduct a brute force search of 4 temperature space and find the
        best fit
        set log to create a log file of the iterations rather that outputting to
        the jupyter/ipython session

        Parameters
        ----------

        initial:  `list`
            the initial trial values (4) for the log10 emission measure


        Keyword Arguments
        -----------------

        indxlimits:  `list`, None
            the range of indices of the density array to search
            if None is specified, the whole range is searched

        verbose:  `bool`
            if True, additional output is sent to the terminal

        log:  `bool`
            if True, a log file is created - 'search1d.raw'

        '''
        self.Nparams = 8. + 1.
        t1 = datetime.now()
        if indxlimits == None:
            if not hasattr(self, 'MinIndex'):
                self.findMinMaxIndices()
            indxlimits = [self.MinIndex, self.MaxIndex]
        print('indxlimits = %5i %5i '%(indxlimits[0], indxlimits[1]))
        nTemp = len(self.EmIndices)
        emfit = []
        idx123 = []
        searchDx123 = []
        chisq = []
        maskedValues = []
        temp1Searched = self.Temperature[indxlimits[0]:indxlimits[1]-2]
        temp2Searched = self.Temperature[indxlimits[0]+1:indxlimits[1]-1]
        temp3Searched = self.Temperature[indxlimits[0]+2:indxlimits[1]]
        temp4Searched = self.Temperature[indxlimits[0]+3:indxlimits[1]+1]

        temp1Size = temp1Searched.size
        temp2Size = temp2Searched.size
        temp3Size = temp3Searched.size
        temp4Size = temp4Searched.size

#        msk  = np.ma.make_mask(np.ones((temp1Searched.size, temp2Searched.size, temp3Searched.size)))
#        msk2 = np.ma.make_mask(np.ones((temp1Searched.size, temp2Searched.size, temp3Searched.size, 3)))
        chisqNd = np.ma.array(np.zeros((temp1Size, temp2Size, temp3Size, temp4Size), 'float64'), mask=True)
        tempNd = np.array(np.zeros((4, temp1Size, temp2Size, temp3Size, temp4Size), 'float64'))
        emNd = np.ma.array(np.zeros((4, temp1Size, temp2Size, temp3Size, temp4Size), 'float64'), mask=True)

        if log:
            logname = 'search4t.raw'
            logfile = open(logname, 'w')
            logfile.write('%s \n'%(self.SpecData['filename']))
        # need to make this look like search2t
        pformat = '%5i %4i %4i %4i %4i %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e'
        logformat = '%5i %4i %4i %4i %4i %4i %4i %4i %4i %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n'
        counter = 0
        for idx1 in range(indxlimits[0], indxlimits[1] - 2, 1):
            for idx2 in range(idx1+1, indxlimits[1] - 1, 1):
                for idx3 in range(idx2+2, indxlimits[1], 1):
                    for idx4 in range(idx3+3, indxlimits[1] + 1, 1):
                        self.emSetIndices([idx1, idx2, idx3, idx4])
                        self.fitNt(initial)
            # need to make this look like search 2t
                        chisq1,  msk1 = self.getChisq()
                        searchDx1 = idx1 - indxlimits[0]
                        searchDx2 = idx2 - indxlimits[0] - 1
                        searchDx3 = idx3 - indxlimits[0] - 2
                        searchDx4 = idx4 - indxlimits[0] - 3
                        searchDx123.append([searchDx1, searchDx2, searchDx3, searchDx4])
                        tempNd[0, searchDx1, searchDx2, searchDx3, searchDx4] = temp1Searched[searchDx1]
                        tempNd[1, searchDx1, searchDx2, searchDx3, searchDx4] = temp2Searched[searchDx2]
                        tempNd[2, searchDx1, searchDx2, searchDx3, searchDx4] = temp3Searched[searchDx3]
                        tempNd[3, searchDx1, searchDx2, searchDx3, searchDx4] = temp4Searched[searchDx4]
                        if not msk1:
                            emfit.append(self.Leastsq['em'])
                            chisq.append(chisq1)
                            idx123.append([idx1, idx2, idx3, idx4])
                            chisqNd.data[searchDx1, searchDx2, searchDx3, searchDx4] = chisq1
                            chisqNd.mask[searchDx1, searchDx2, searchDx3, searchDx4] = False
                            for it in range(nTemp):
                                emNd.data[it, searchDx1, searchDx2, searchDx3, searchDx4] = self.Leastsq['em'][it]
                                emNd.mask[it, searchDx1, searchDx2, searchDx3, searchDx4] = False

                            if verbose:
                                print(pformat%(counter, idx1, idx2, idx3, idx4, self.Temperature[idx1], self.Leastsq['em'][0], self.Temperature[idx2], self.Leastsq['em'][1], self.Temperature[idx3], self.Leastsq['em'][2], self.Temperature[idx4], self.Leastsq['em'][3], chisq[-1]))
                            elif log:
                                logfile.write(logformat%(counter, idx1, idx2, idx3, idx4, searchDx1, searchDx2, searchDx3, searchDx4, self.Temperature[idx1], self.Leastsq['em'][0], self.Temperature[idx2], self.Leastsq['em'][1], self.Temperature[idx3], self.Leastsq['em'][2], self.Temperature[idx4], self.Leastsq['em'][3], chisq[-1]))

                        else:
                            emfit.append(self.Leastsq['em'])
                            idx123.append([idx1, idx2, idx3, idx4])
                            chisq.append(chisq1)
                            maskedValues.append([idx1, idx2, idx3, idx4])
                            chisqNd.data[searchDx1, searchDx2, searchDx3, searchDx4] = chisq1
                            chisqNd.mask[searchDx1, searchDx2, searchDx3, searchDx4] = True
                            for it in range(nTemp):
                                emNd.data[it, searchDx1, searchDx2, searchDx3, searchDx4] = self.Leastsq['em'][it]
                                emNd.mask[it, searchDx1, searchDx2, searchDx3, searchDx4] = True

                            if verbose:
                                print('msk at %5i %4i %4i %4i %4i %10.2e'%(counter, idx1, idx2, idx3, idx4, chisq[-1]))
                            elif log:
                                logfile.write('msk at %5i %4i %4i %4i %4i %10.2e \n'%(counter, idx1, idx2, idx3, idx4, chisq[-1]))
                        counter += 1
        if log:
            logfile.close()
        #
        self.SearchData = {'temperature':self.Temperature,'temp1Searched':temp1Searched, 'temp2Searched':temp2Searched,
        'temp3Searched':temp3Searched, 'temp4Searched':temp4Searched, 'emfit':emfit, 'idx':idx123, 'searchDx123':searchDx123, 'chisq':chisq, 'chisqNd':chisqNd,
        'minchiseq':min(chisq), 'tempNd':tempNd, 'emNd':emNd, 'message':'this contains all the data for the 3tExpSpace', 'nparams':self.Nparams}
        #
        #
        #  set things to best fit
        #
        print('min chisq = %10.3e at %i'%(min(chisq), np.argmin(chisq)))

        gdx=[i for i,ch in enumerate(chisq) if ch == min(chisq)]
        if verbose:
#        print(' gdx = ', gdx)
            print(' len, min of chisq %5i %10.3e'%(len(chisq), min(chisq)))
            if len(gdx) == 1:
                print('single value gdx = %i '%(gdx[0]))  #, idx12[gdx[0]][0], idx12[gdx[0]][1])
#                temperature1 = eis.SearchData['temp1Searched']
            else:
                pstring = 'mult values of gdx = '
                for one in gdx:
                    pstring += '  %10i'%(one)
        # 90% confidence for 6 parameters
#        chisq1sig = min(chisq)*np.ones_like(chisqNd, 'float64') + 13.36
#        good1sig = chisqNd < min(chisq) + 13.36

        self.emSetIndices(idx123[gdx[0]])
        self.emSet(emfit[gdx[0]])
        self.predict()
        em = 10.**emfit[gdx[0]]
        emfit1 = emfit[gdx[0]]
        adx = idx123[gdx[0]]
        minChisq = chisq[gdx[0]]
        reducedChisq = minChisq/float(self.Nobs - 2*nTemp)
        print('reduced Chisq = %12.4e'%(reducedChisq))
        print(' idx1 %5i idx2 %5i idx3 %5i idx4 %5i '%(adx[0], adx[1], adx[2], adx[3]))
        print(' Temp =  %10.2e %10.2e %10.2e %10.2e'%(self.Temperature[adx[0]], self.Temperature[adx[1]], self.Temperature[adx[2]], self.Temperature[adx[3]]))
        print(' EM =    %10.2e %10.2e %10.2e %10.2e'%(em[0], em[1], em[2], em[3]))
        print(' EmFit = %10.3f %10.3f %10.3f %10.3f'%(emfit1[0], emfit1[1], emfit1[2], emfit1[3]))
        #
        t2 = datetime.now()
        dt = t2 - t1
        print(' elapsed seconds = %12.3f'%(dt.seconds))

        self.SearchData['best'] = {'em':em,'emfit':emfit[gdx[0]], 'chisq':minChisq, 'reducedChisq': reducedChisq, 'idx':idx123[gdx[0]], 'gdx':gdx[0],
        'temperature':self.Temperature[idx123[gdx[0]]], 'density':self.Density}


    def loadMatch(self, filename):
        """ to open a pickle file, return the match data and make it an attribute
        """
        with open(filename, 'rb') as inpt:
            matchDict = pickle.load(inpt)
        self.Match = matchDict['match']
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

        if 'EmLog' in matchDict.keys():
            self.EmLog = matchDict['EmLog']
        else:
            print('Em not in matchDict')

        if 'NT' in matchDict.keys():
            self.NT = matchDict['NT']
        else:
            print('NT not in matchDict')

        if 'XUVTOP' in matchDict.keys():
            self.XUVTOP = matchDict['XUVTOP']
        if 'chiantiVersion' in matchDict.keys():
            self.ChiantiVersion = matchDict['chiantiVersion']
        self.MatchName = filename
        #
        #-----------------------------------------------------
        #
    def saveMatch(self, filename):
        """to save the attribute Match to a pickle file so that it can be reloaded later


        Keyword Arguments
        -----------------

        filename:  `str`
            the filename where the text should be output

        """
        matchDict={'match':self.Match, 'Temperature':self.Temperature, 'EDensity':self.EDensity, 'Ndens':self.Ndens,
            'Ntemp':self.Ntemp, 'NTempDens':self.NTempDens, 'MinAbund':self.MinAbund}
        if hasattr(self, 'EmIndices'):
            matchDict['EmIndices'] = self.EmIndices

        if hasattr(self, 'Em'):
            matchDict['Em'] = self.Em

        if hasattr(self, 'EmLog'):
            matchDict['EmLog'] = self.EmLog

        if hasattr(self, 'Nparams'):
            matchDict['Nparams'] = self.Nparams

        if hasattr(self, 'NT'):
            matchDict['NT'] = self.NT
        else:
            print('NT not available')

        if hasattr(self, 'WghtFactor'):
            matchDict['WghtFactor'] = self.WghtFactor
        else:
            print('wghtfactor not available')

        if hasattr(self, 'AbundanceName'):
            matchDict['AbundanceName'] = self.AbundanceName
        else:
            print('AbundanceName not available')

        if 'XUVTOP' in self.SpecData.keys():
            matchDict['XUVTOP'] = self.SpecData['XUVTOP']

        if 'chiantiVersion' in self.SpecData.keys():
            matchDict['chiantiVersion'] = self.SpecData['chiantiVersion']

        with open(filename, 'wb') as outpt:
            pickle.dump(matchDict, outpt)

    def loadSearchData(self, filename):
        """ to load the pickled search data as an attribute self.SearchData


        Keyword Arguments
        -----------------

        filename:  `str`
            the filename of the pickle file where the search data has been created


        """
        with open(filename, 'rb') as inpt:
            self.SearchData = pickle.load(inpt)

    def saveSearchData(self, filename):
        """to save the attribute SearchData to a pickle file


        Keyword Arguments
        -----------------

        filename:  `str`
            the filename of the pickle file where the search data is to be created


        """
        with open(filename, 'wb') as outpt:
            pickle.dump(self.SearchData, outpt)
