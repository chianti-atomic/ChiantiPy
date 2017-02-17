"""
Base class used in several ChiantiPy objects
"""

from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

import ChiantiPy.tools.filters as chfilters
import ChiantiPy.tools.util as util
import ChiantiPy.tools.io as chio
import ChiantiPy.tools.data as chdata
import ChiantiPy.tools.constants as const


class specTrails(object):
    """
    a collection of methods for use in spectrum calculations
    """


    def __init__(self, temperature, density):
        self.Temperature = temperature
        self.EDensity = density
        self.AbundanceName = chdata.Defaults['abundfile']
        self.AbundAll = chdata.Abundance[self.AbundanceName]['abundance']
        #
        # ---------------------------------------------------------------------------
        #
    def convolve(self, wavelength=0, filter=(chfilters.gaussianR, 1000.), label=0, verbose=0):
        '''
        the first application of spectrum calculates the line intensities within the specified wavelength range and for set of ions specified

        wavelength will not be used if applied to 'spectrum' objects

        wavelength IS need for 'bunch' objects - in this case, the wavelength should not extend beyond the limits of the
        wvlRange used for the 'bunch' calculation

        '''
        if not hasattr(self, 'IonInstances'):
            print(' must set keepIons=1 in order to keep self.IonInstances')
            return
        #
        if type(label)!= type(0):
            if type(label) != str:
                print(' label must either be zero or a string')
                return
        #
        t1 = datetime.now()
        #:
        if hasattr(self, 'Wavelength'):
                nWvl = len(self.Wavelength)
                wavelength = self.Wavelength
        elif type(wavelength) == int:
            print(' a wavelength array must be given')
            return
        else:
            self.Wavelength = wavelength
            nWvl = len(wavelength)
        lineSpectrum = np.zeros((self.NTempDen, nWvl), 'float64').squeeze()
        for akey in sorted(self.IonInstances.keys()):
            if verbose:
                print( ' trying ion = %s'%(akey))
#            thisIon = self.IonInstances[akey]
            if not 'errorMessage' in sorted(self.IonInstances[akey].Intensity.keys()):
                if verbose:
                    print(' doing convolve on ion %s '%(akey))
                self.IonInstances[akey].spectrum(wavelength, filter)
#                lineSpectrum = np.add(lineSpectrum, self.IonInstances[akey].Spectrum['intensity'])
                if 'errorMessage' in sorted(self.IonInstances[akey].Spectrum.keys()):
                    print(self.IonInstances[akey].Spectrum['errorMessage'])
                else:
                    lineSpectrum += self.IonInstances[akey].Spectrum['intensity']
#                if self.NTempDen == 1:
#                    lineSpectrum += thisIon.Spectrum['intensity']
#                else:
#                    for iTempDen in range(self.NTempDen):
#                        lineSpectrum[iTempDen] += thisIon.Spectrum['intensity'][iTempDen]
            else:
                if 'errorMessage' in sorted(self.IonInstances[akey].Intensity.keys()):
                    print(self.IonInstances[akey].Intensity['errorMessage'])

        self.LineSpectrum = {'wavelength':wavelength, 'intensity':lineSpectrum.squeeze()}
        #
        total = self.LineSpectrum['intensity']
        #
        # the following is required in order to be applied to both a 'spectrum' and a 'bunch' object
        #
        if hasattr(self, 'FreeFree'):
            total += self.FreeFree['intensity']
        if hasattr(self, 'FreeBound'):
            total += self.FreeBound['intensity']
        if hasattr(self, 'TwoPhoton'):
            total += self.TwoPhoton['intensity']
        self.Total = total
        #
        #
        if self.NTempDen == 1:
            integrated = total
        else:
            integrated = total.sum(axis=0)
        #
        t2 = datetime.now()
        dt = t2 - t1
        print(' elapsed seconds = %12.3e'%(dt.seconds))
        #
        if type(label) == type(''):
            if hasattr(self, 'Spectrum'):
                self.Spectrum[label] = {'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':self.Em,  'Abundance':self.AbundanceName}
            else:
                self.Spectrum = {label:{'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':self.Em,  'Abundance':self.AbundanceName}}
        else:
            self.Spectrum ={'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'Abundance':self.AbundanceName}
        return
        #
        # ---------------------------------------------------------------------------
        #
    def ionGate(self, elementList = 0, ionList = 0, minAbund=0, doContinuum=1, doWvlTest=1, verbose=0):
        '''
        creates a list of ions for free-free, free-bound, and line intensity calculations
        if doing the radiative losses, accept all wavelength -> doWvlTest=0
        the list is a dictionary self.Todo
        '''
        #
        masterlist = chdata.MasterList
        abundAll = self.AbundAll
        #
        nonzed = abundAll > 0.
        minAbundAll = abundAll[nonzed].min()
        if minAbund:
            if minAbund < minAbundAll:
                minAbund = minAbundAll
        ionInfo = chio.masterListInfo()
        #
        if hasattr(self, 'Wavelength'):
            wvlRange = [self.Wavelength.min(), self.Wavelength.max()]
        elif hasattr(self, 'WvlRange'):
            wvlRange = self.WvlRange
        else:
            print(' need a wavelength range in ionGate ')
        #
        temperature = self.Temperature
        #
        # use the ionList but make sure the ions are in the database
        self.Todo = {}
        #
        if minAbund:
            if doContinuum:
                for iz in range(1, 31):
                    abundance = chdata.Abundance[self.AbundanceName]['abundance'][iz-1]
                    if abundance >= minAbund:
                        if verbose:
                            print(' %5i %5s abundance = %10.2e '%(iz, const.El[iz-1],  abundance))
                        #
                        for ionstage in range(1, iz+2):
                            ionS = util.zion2name(iz, ionstage)
                            masterListTest = ionS in masterlist
                            masterListInfoTest = ionS in sorted(ionInfo.keys())
                            if masterListTest or masterListInfoTest:
                                if doWvlTest:
                                    wvlTestMin = wvlRange[0] <= ionInfo[ionS]['wmax']
                                    wvlTestMax = wvlRange[1] >= ionInfo[ionS]['wmin']
                                else:
                                    wvlTestMin = 1
                                    wvlTestMax = 1
                                ioneqTest = (temperature.max() >= ionInfo[ionS]['tmin']) and (temperature.min() <= ionInfo[ionS]['tmax'])
                            # construct similar test for the dielectronic files
                            ionstageTest = ionstage > 1
                            if ionstageTest and ioneqTest and doContinuum:
                                # ionS is the target ion, cannot be the neutral for the continuum
                                if verbose:
                                    print(' setting up continuum calculation for %s  '%(ionS))
                                if ionS in sorted(self.Todo.keys()):
                                    self.Todo[ionS] += '_ff_fb'
                                else:
                                    self.Todo[ionS] = 'ff_fb'
                                if verbose:
                                    print(' for ion %s do : %s'%(ionS, self.Todo[ionS]))
        #
        if elementList:
            for i,  one in enumerate(elementList):
                elementList[i] = one.lower()
            for one in masterlist:
                stuff = util.convertName(one)
                bare = stuff['Z'] == stuff['Ion']
                if stuff['Element'] in  elementList:
                    self.Todo[one] = 'line'
                    if doContinuum and not stuff['Dielectronic']:
                        self.Todo[one]+= '_ff'
                        if not bare:
                            self.Todo[one] += '_fb'
        if ionList:
            for one in ionList:
                stuff = util.convertName(one)
                bare = stuff['Z'] == stuff['Ion']
                if masterlist.count(one):
                    self.Todo[one] = 'line'
                    if doContinuum and not stuff['Dielectronic']:
                        self.Todo[one]+= '_ff'
                        if not bare:
                            self.Todo[one] += '_fb'
                else:
                    if verbose:
                        pstring = ' %s not in CHIANTI database'%(one)
                        print(pstring)
        #
        #
        #
        if minAbund:
            for iz in range(1, 31):
                abundance = chdata.Abundance[self.AbundanceName]['abundance'][iz-1]
                if abundance >= minAbund:
                    if verbose:
                        print(' %5i %5s abundance = %10.2e '%(iz, const.El[iz-1],  abundance))
                    #
                    for ionstage in range(1, iz+2):
                        ionS = util.zion2name(iz, ionstage)
                        masterListTest = ionS in masterlist
                        masterListInfoTest = ionS in sorted(ionInfo.keys())
                        if masterListTest or masterListInfoTest:
                            if masterListTest or masterListInfoTest:
                                if doWvlTest:
                                    wvlTestMin = wvlRange[0] <= ionInfo[ionS]['wmax']
                                    wvlTestMax = wvlRange[1] >= ionInfo[ionS]['wmin']
                                else:
                                    wvlTestMin = 1
                                    wvlTestMax = 1
                            ioneqTest = (temperature.max() >= ionInfo[ionS]['tmin']) and (temperature.min() <= ionInfo[ionS]['tmax'])
                        # construct similar test for the dielectronic files
                        ionSd = util.zion2name(iz, ionstage, dielectronic=1)
                        masterListTestD = ionSd in masterlist
                        masterListInfoTestD = ionSd in sorted(ionInfo.keys())
                        if masterListTestD or masterListInfoTestD:
                            if doWvlTest:
                                wvlTestMinD = wvlRange[0] <= ionInfo[ionSd]['wmax']
                                wvlTestMaxD = wvlRange[1] >= ionInfo[ionSd]['wmin']
                            else:
                                wvlTestMinD = 1
                                wvlTestMaxD = 1
                            ioneqTestD = (temperature.max() >= ionInfo[ionSd]['tmin']) and (temperature.min() <=ionInfo[ionSd]['tmax'])
                            #
                        if masterListTest and wvlTestMin and wvlTestMax and ioneqTest:
                            #if verbose:
                                #print(' setting up spectrum calculation for  %s'%(ionS))
                            if ionS in sorted(self.Todo.keys()):
                                self.Todo[ionS] += '_line'
                            else:
                                self.Todo[ionS] = 'line'
                        # get dielectronic lines
                            if verbose:
                                print(' for ion %s do : %s'%(ionS, self.Todo[ionS]))
                        if masterListTestD and wvlTestMinD and wvlTestMaxD and ioneqTestD:
                            #if verbose:
                                #print(' setting up  spectrum calculation for  %s '%(ionSd))
                            if ionSd in sorted(self.Todo.keys()):
                                self.Todo[ionSd] += '_line'
                            else:
                                self.Todo[ionSd] = 'line'
                            if verbose:
                                print(' for ion %s do : %s'%(ionSd, self.Todo[ionSd]))
        if len(self.Todo.keys()) == 0:
            print(' no elements have been selected')
            print(' it is necessary to provide an elementList, an ionList, or set minAbund > 0.')
        return
    #
    # -------------------------------------------------------------------------
    #
    def spectrumPlot(self, index=-1, integrated=0, saveFile=0, linLog = 'lin'):
        '''
        to plot the spectrum as a function of wavelength
        '''
        plt.figure()
        mask = self.Em > 1.
        if mask.sum() == 0:
            ylabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1} \AA^{-1}$ ($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'
        else:
            ylabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1} \AA^{-1}$'
        #
        xlabel = 'Wavelength ('+self.Defaults['wavelength'] +')'
        #
#        ymin = 10.**(np.log10(emiss.min()).round(0))
        #
        plt.ion()
        #
        if integrated:
            if 'wavelength' in sorted(self.Spectrum.keys()):
                plt.plot(self.Spectrum['wavelength'], self.Spectrum['integrated'])
            elif 'wvl' in sorted(self.Spectrum.keys()):
                plt.plot(self.Spectrum['wvl'], self.Spectrum['integrated'])
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title('integrated spectrum')
        else:
            nTempDen = self.NTempDen
            if nTempDen == 1:
            #
                if 'wavelength' in sorted(self.Spectrum.keys()):
                    plt.plot(self.Spectrum['wavelength'], self.Spectrum['intensity'])
                elif 'wvl' in sorted(self.Spectrum.keys()):
                    plt.plot(self.Spectrum['wvl'], self.Spectrum['intensity'])
                    plt.title(' Temperature = %10.2e K'%(self.Temperature))
            else:
                if index < 0:
                    index = nTempDen/2
                if 'wavelength' in sorted(self.Spectrum.keys()):
                    plt.plot(self.Spectrum['wavelength'], self.Spectrum['intensity'][index])
                elif 'wvl' in sorted(self.Spectrum.keys()):
                    plt.plot(self.Spectrum['wvl'], self.Spectrum['intensity'][index])
                    plt.title(' Temperature = %10.2e K for index = %3i'%(self.Temperature[index], index))
                plt.xlabel(xlabel)
                plt.ylabel(ylabel)
        if saveFile:
            plt.savefig(saveFile)
    #
    # -------------------------------------------------------------------------
    #
    def lineSpectrumPlot(self, index = 0, integrated=0, saveFile=0, linLog = 'lin'):
        '''
        to plot the line spectrum as a function of wavelength
        '''
        #
        #
        plt.figure()
        mask = self.Em > 1.
        if mask.sum() == 0:
            ylabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1} \AA^{-1}$ ($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'
        else:
            ylabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1} \AA^{-1}$'

        #
        xlabel = 'Wavelength ('+self.Defaults['wavelength'].capitalize() +')'
        #
#        ymin = 10.**(np.log10(emiss.min()).round(0))
        #
        plt.ion()
        #
        if integrated:
            plt.plot(self.Spectrum['wavelength'], self.Spectrum['intensity'].sum(axis=0))
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title('integrated spectrum')
        else:
            nTempDen = self.NTempDen
            if nTempDen == 1:
            #
                plt.plot(self.LineSpectrum['wavelength'], self.LineSpectrum['intensity'])
                plt.title(' Temperature = %10.2e K '%(self.Temperature))
            else:
                plt.plot(self.LineSpectrum['wavelength'], self.LineSpectrum['intensity'][index])
                plt.title(' Temperature = %10.2e K for index = %3i'%(self.Temperature[index], index))
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
        if saveFile:
            plt.savefig(saveFile)
    #
    # -------------------------------------------------------------------------
