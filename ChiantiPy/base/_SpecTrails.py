"""
Base class used in several ChiantiPy objects
"""

from datetime import datetime
import pickle

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
    def convolve(self, wavelength=0, filter=(chfilters.gaussianR, 1000.), label=0, verbose=False):
        '''
        the first application of spectrum calculates the line intensities within the specified wavelength range and for set of ions specified

        wavelength will not be used if applied to 'spectrum' objects

        wavelength IS need for 'bunch' objects - in this case, the wavelength should not extend beyond the limits of the
        wvlRange used for the 'bunch' calculation

        Keyword Arguments
        -----------------

        wavelength:  'int', `list`
            if an `int`, the attribute 'Wavelength' is looked for
            otherwise, wavelength is used

        filter: `tuple`
            first elements if one of the ChiantiPy.tools.filters object
            second element is the width appropriate to the filter

        label:  `str`
            if set, creates a Spectrum[label] attribute

        verbose: `bool`
            if True, prints info to the terminal

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
        lineSpectrum = np.zeros((self.NTempDens, nWvl), np.float64).squeeze()
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
#                if self.NTempDens == 1:
#                    lineSpectrum += thisIon.Spectrum['intensity']
#                else:
#                    for iTempDen in range(self.NTempDens):
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
        if self.NTempDens == 1:
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
    def ionGate(self, elementList = None, ionList = None, minAbund=None, doLines=1, doContinuum=1, doWvlTest=1, doIoneqTest=1, includeDiel=False,  verbose=0):
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
        #
        todo = {}
        zTodo = []
        #
        if elementList:
            elementListLower = []
            for i,  element in enumerate(elementList):
                elementListLower.append(element.lower())
        #
        if elementList:
            if verbose:
                print(' \n in elementList \n')
            for i,  element in enumerate(elementListLower):
                if verbose:
                    print('el = %s'%(element))
                z = util.el2z(element)
                zTodo.append(z)
                for stage in range(1, z+1):
                    ion = util.zion2name(z,  stage)
                    if ion in masterlist:
                        nameDict = util.convertName(ion)
                        if doLines:
                                todo[ion] = 'line'
                        if doContinuum and not nameDict['Dielectronic']:
                            # bare ions will be added later
                            if stage != 1:
                                if ion not in todo.keys():
                                    todo[ion] = 'ff'
                                else:
                                    todo[ion] += '_ff_fb'
                        if verbose:
                            print(' %s  %s '%(ion,  todo[ion]))
        if verbose:
            print(' \n after elementList \n')
            for anion in todo.keys():
                print(' %s  %s'%(anion,  todo[anion]))

        if ionList:
            if verbose:
                print(' \n in ionList \n')
            for one in ionList:
                nameDict = util.convertName(one)
                zTodo.append(nameDict['Z'])
                if masterlist.count(one):
                    if doLines:
                        todo[one] = 'line'
                    if verbose:
                        print(' %s in the CHIANTI database'%(one))
                else:
                    if verbose:
                        pstring = ' %s not in CHIANTI database'%(one)
                        print(pstring)
                if doContinuum and not nameDict['Dielectronic'] and nameDict['Ion'] < 31:
                    if one not in todo.keys():
                        todo[one] = 'ff_fb'
                    else:
                        todo[one] += '_ff_fb'

        # add bare nuclei
        if doContinuum:
            zSet = set(zTodo)
            for z in zSet:
                gname = util.zion2name(z, z+1)
                todo[gname] = 'ff_fb'


        if verbose:
            print(' \n after ionList \n')
            for anion in todo.keys():
                print(' %s  %s'%(anion,  todo[anion]))


        if doIoneqTest:
            if verbose:
                print(' \n in ioneq test \n')
            toPop = []
            for ionS in todo.keys():
                ioneqTest = (temperature.max() >= ionInfo[ionS]['tmin']) and (temperature.min() <= ionInfo[ionS]['tmax'])
                if not ioneqTest:
                    toPop.append(ionS)
                else:
                    if verbose:
                        print('passes ioneqtest  %s'%( ionS))
            for badion in toPop:
                if verbose:
                    print('fails ioneqtest %s'%(badion))
                todo.pop(badion)

        if verbose:
            print(' \n after ioneq test \n')
            for anion in todo.keys():
                print(' %s  %s'%(anion,  todo[anion]))

        if doWvlTest:
            if verbose:
                print(' \n doing wvltest \n')
            toPop = []
            for ionS in todo:
                # bare ions have wmin=0 and wmax=1.e+30
                if doWvlTest:
                    wvlTestMin = wvlRange[0] <= ionInfo[ionS]['wmax']
                    wvlTestMax = wvlRange[1] >= ionInfo[ionS]['wmin']
                else:
                    wvlTestMin = 1
                    wvlTestMax = 1
#                if verbose:
#                    print(' %s  %8.2f  %8.2f %8.2f %8.2f'%(ionS, ionInfo[ionS]['wmin'], ionInfo[ionS]['wmax'], wvlRange[0], wvlRange[1]))
#                    print(' %s wvlTestMin  %s  wvlTestMax %s'%(ionS, wvlTestMin, wvlTestMax))
                if wvlTestMin and wvlTestMax:
                    if verbose:
                        print(' %s passed wvl test'%(ionS))
                else:
                    toPop.append(ionS)
            for badion in toPop:
                if verbose:
                    print('fails wvltest %s'%(badion))
            for badion in toPop:
                todo[badion] = todo[badion].replace('line', '')

        if verbose:
            print(' \n after wvl test \n')
            for anion in todo.keys():
                print(' %s  %s'%(anion,  todo[anion]))
        #
        #
        #  the relative H abundance is 1.0, the rest are all smaller
        if minAbund is not None and type(minAbund) is float:
            if verbose:
                print(' \n doing minAbund test \n')
            for iz in range(1, 31):
#                abundance = chdata.Abundance[self.AbundanceName]['abundance'][iz-1]
                abundance = self.Abundance[iz-1]
                if abundance >= minAbund:
                    if verbose:
                        print(' %5i %5s abundance = %10.2e '%(iz, const.El[iz-1],  abundance))
                    #
                    for ionstage in range(1, iz+1):
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
                            if temperature.size == 1:
                                ioneqTest = (temperature >= ionInfo[ionS]['tmin']) and (temperature <= ionInfo[ionS]['tmax'])
                            else:
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
                        if masterListTest and wvlTestMin and wvlTestMax and ioneqTest and doLines:
                            if verbose:
                                print(' %s passed mList, wvlTest, ioneqTest'%(ionS))
                            if ionS in sorted(todo.keys()):
                                todo[ionS] += '_line'
                            else:
                                todo[ionS] = 'line'
                        # get dielectronic lines
                            if verbose:
                                print(' for ion %s do : %s'%(ionS, todo[ionS]))
                        if masterListTestD and wvlTestMinD and wvlTestMaxD and ioneqTestD and doLines:
                            if ionSd in sorted(todo.keys()):
                                todo[ionSd] += '_line'
                            else:
                                todo[ionSd] = 'line'
                            if verbose:
                                print(' for ion %s do : %s'%(ionSd, todo[ionSd]))
        #
                    #
                    for ionstage in range(2, iz+2):
                        ionS = util.zion2name(iz, ionstage)
                        if ionS in ionInfo.keys():
                            ioneqTest = (temperature.max() >= ionInfo[ionS]['tmin']) and (temperature.min() <= ionInfo[ionS]['tmax'])
                        else:
                            ioneqTest = 1
                        # construct similar test for the dielectronic files
                        if ioneqTest and doContinuum:
                            # ionS is the target ion, cannot be the neutral for the continuum
#                            if verbose:
#                                print(' setting up continuum calculation for %s  '%(ionS))
                            if ionstage == iz+1 and doContinuum:
                                # this is the bare ion
                                if ionS in sorted(todo.keys()):
                                    todo[ionS] += '_ff_fb'
                                else:
                                    todo[ionS] = 'ff_fb'
                                if verbose:
                                    print(' for ion %s do : %s'%(ionS, todo[ionS]))
                            elif doContinuum:
                                if ionS in sorted(todo.keys()):
                                    todo[ionS] += '_ff_fb'
                                else:
                                    todo[ionS] = 'ff_fb'
                                if verbose:
                                    print(' for ion %s do : %s'%(ionS, todo[ionS]))
        if verbose:
            print(' finished minabund test \n')

        if verbose:
            print(' \n after minabund test \n')
            for anion in todo.keys():
                print(' %s  %s'%(anion,  todo[anion]))

        #  remove dupicates
#        todoSet = set(todo)
#        self.Todo = list(todoSet)
        dielList = []
        if not includeDiel:
            for akey in todo:
                if 'd' in akey[-1]:
                    dielList.append(akey)
                    if verbose:
                        print(' added dielectronic ion to diellist%s'%(akey))
        newTodo = {}
        for akey in todo:
            if akey not in dielList:
                newTodo[akey] = todo[akey]
        self.Todo = newTodo
        if len(self.Todo.keys()) == 0:
            print(' no elements have been selected')
            print(' it is necessary to provide an elementList, an ionList, or set minAbund > 0.')
        return
    #
    # -------------------------------------------------------------------------
    #
    def saveData(self, filename):
        """

        :param filename: filename where the pickle file of the saved data will be stored
        :type filename: str

        following running a ch.spectrum calculation, save the calculation as a dictionary to a pickle file
        """
        data = {'temperature':self.Temperature, 'eDensity':self.EDensity,
            'em':self.Em, 'abundanceName':self.AbundanceName, 'abundAll':self.AbundAll,
            'ionsCalculated':self.IonsCalculated,
            'defaults':self.Defaults, 'intensity':self.Intensity,
            'nTemp':self.Ntemp, 'nDens':self.Ndens, 'nTempDens':self.NTempDens}

#          'elementList':self.elementList, 'ionList':self.ionList,
#            'minAbund':self.minAbund, 'keepIons':self.keepIons, 'em':self.em, 'abundance':self.abundance,
#            'allLines':self.allLines}

        if hasattr(self, 'Spectrum'):
            data['spectrum'] = self.Spectrum
        if hasattr(self, 'Wvl'):
            data['wvl'] = self.Wvl
        if hasattr(self, 'Wavelength'):
            data['wavelength'] = self.Wavelength
        if hasattr(self, 'WvlRange'):
            data['wvlRange'] = self.WvlRange
        if hasattr(self, 'PDensity'):
            data['pDensity'] = self.PDensity
        if hasattr(self, 'IonInstances'):
            data['ionInstances'] = self.IonInstances
        if hasattr(self, 'Xlabel'):
            data['xlabel'] = self.Xlabel
        if hasattr(self, 'Ylabel'):
            data['ylabel'] = self.Ylabel
        with open(filename, 'wb') as outpt:
            pickle.dump(data, outpt)

    def spectrumPlot(self, index=-1, integrated=False, saveFile=False, linLog = 'lin'):
        '''
        to plot the spectrum as a function of wavelength

        Keyword Arguments
        -----------------

        index:  `int`
            selects the temperature of the calculated spectrum

        integrated:  `bool`
            if True, plots the integrated/summed spectrum

        saveFile:  `bool`, `str`
            if set, saves the plot to 'saveFile'

        linLog:  `str`
            should be either 'lin' for linear, or 'log' for logarithmic base10

        '''
        fs = 14
        plt.figure()
        mask = self.Em > 1.
        if mask.sum() == 0:
            self.Ylabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1} \AA^{-1}$ ($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'
        else:
            self.Ylabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1} \AA^{-1}$'
        #
        self.Xlabel = 'Wavelength ('+self.Defaults['wavelength'] +')'
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
            plt.xlabel(self.Xlabel,  fontsize=fs)
            plt.ylabel(self.Ylabel,  fontsize=fs)
            plt.title('integrated spectrum',  fontsize=fs)
        else:
            nTempDens = self.NTempDens
            if nTempDens == 1:
            #
                if 'wavelength' in sorted(self.Spectrum.keys()):
                    plt.plot(self.Spectrum['wavelength'], self.Spectrum['intensity'])
                elif 'wvl' in sorted(self.Spectrum.keys()):
                    plt.plot(self.Spectrum['wvl'], self.Spectrum['intensity'])
                    plt.title(' Temperature = %10.2e K'%(self.Temperature),  fontsize=fs)
            else:
                if index < 0:
                    index = nTempDens/2
                if 'wavelength' in sorted(self.Spectrum.keys()):
                    plt.plot(self.Spectrum['wavelength'], self.Spectrum['intensity'][index])
                elif 'wvl' in sorted(self.Spectrum.keys()):
                    plt.plot(self.Spectrum['wvl'], self.Spectrum['intensity'][index])
                    plt.title(' Temperature = %10.2e K for index = %3i'%(self.Temperature[index], index),  fontsize=fs)
                plt.xlabel(self.Xlabel,  fontsize=fs)
                plt.ylabel(self.Ylabel,  fontsize=fs)
        ylim = plt.ylim()
        plt.ylim([0., ylim[1]])
        plt.tight_layout()
        if saveFile:
            plt.savefig(saveFile)
    #
    # -------------------------------------------------------------------------
    #
    def lineSpectrumPlot(self, index = 0, integrated=False, saveFile=False, linLog = 'lin'):
        '''
        to plot the line spectrum as a function of wavelength

        Keyword Arguments
        -----------------

        index:  `int`
            selects the temperature of the calculated line spectrum

        integrated:  `bool`
            if True, plots the integrated/summed line spectrum

        saveFile:  `bool`, `str`
            if set, saves the plot to 'saveFile'

        linLog:  `str`
            should be either 'lin' for linear, or 'log' for logarithmic base10

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
            nTempDens = self.NTempDens
            if nTempDens == 1:
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
