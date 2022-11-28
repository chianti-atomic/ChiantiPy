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


#    def __init__(self, temperature, density):
#        self.Temperature = temperature
#        self.EDensity = density
#        self.AbundanceName = chdata.Defaults['abundfile']
#        self.AbundAll = chdata.Abundance[self.AbundanceName]['abundance']
#        #
        # ---------------------------------------------------------------------------
        #
    def ionGate(self, elementList = None, ionList = None, minAbund=None, doLines=1, doContinuum=1,
        doWvlTest=1, doIoneqTest=1, includeDiel=False,  verbose=0):
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
                minAbund = float(minAbundAll)
        ionInfo = chio.masterListInfo()
        #
        if doWvlTest:
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

                if self.Defaults['wavelength'] == 'angstrom':
                    wmin = ionInfo[ionS]['wmin']
                    wmax = ionInfo[ionS]['wmax']
                elif self.Defaults['wavelength'] == 'nm':
                    wmin = ionInfo[ionS]['wmin']/10.
                    wmax = ionInfo[ionS]['wmax']/10.
                elif self.Defaults['wavelength'] == 'ev':
                    wmin = const.ev2Ang/ionInfo[ionS]['wmax']
                    wmax = const.ev2Ang/ionInfo[ionS]['wmin']
                elif self.Defaults['wavelength'] == 'kev':
                    wmin = const.kev2Ang/ionInfo[ionS]['wmax']
                    wmax = const.kev2Ang/ionInfo[ionS]['wmin']

                if verbose:
                    print('%s wmin:  %10.2e  wmax:  %10.2e'%(ionS,  wmin,  wmax))

                if doWvlTest:
                    wvlTestMin = wvlRange[0] <= wmax
#                    print(' wvlRange[0] %10.2e <= wmax %10.2e'%(wvlRange[0],  wmax))
                    wvlTestMax = wvlRange[1] >= wmin
#                    print(' wvlRange[1] %10.2e >= wmin %10.2e'%(wvlRange[1],  wmin))
                else:
                    wvlTestMin = 1
                    wvlTestMax = 1
                if verbose:
                    print(' %s  %8.2f  %8.2f %8.2f %8.2f'%(ionS, ionInfo[ionS]['wmin'], ionInfo[ionS]['wmax'], wvlRange[0], wvlRange[1]))
                    print(' %s wvlTestMin  %s  wvlTestMax %s'%(ionS, wvlTestMin, wvlTestMax))
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
    def restoreData(self,  filename):
        """

        :param filename: filename where the pickle file of the saved data to be loaded
        :type filename: str

        """
        with open(filename, 'rb') as inpt:
            data = pickle.load(inpt)

        for akey in data.keys():
            setattr(self, akey,  data[akey])
        return

    #
    # -------------------------------------------------------------------------
    #
    def saveData(self, filename, verbose=False):
        """

        :param filename: filename where the pickle file of the saved data will be stored
        :type filename: str

        following running a multi-ion calculation (bunch, spectrum, mspectrum, ipymspectrum, radloss,
        save the calculation as a dictionary to a pickle file
        """

        t0 = datetime.now()
        data = {'Filename':filename, 'Date':t0.strftime('%Y %B %d %H%M'),
            'ClassName':self.__class__.__name__}
        for aname in self.__dict__.keys():
            if hasattr(self, aname):
                data[aname] = getattr(self, aname)

        with open(filename, 'wb') as outpt:
            pickle.dump(data, outpt)

        if verbose:
            for aname  in data:
                print(' saving attribute:  %s'%(aname))



    def spectrumPlot(self,  wvlRange=None, index=None, integrated=False, saveFile=False, linLog = 'lin',
        doLabel=True, lw=1, doTitle=True, top=10):
        '''
        to plot the spectrum as a function of wavelength

        Keyword Arguments
        -----------------

        wvlRange:  2 element tuple, list or array determines the wavelength range to plot

        index:  `int`
            selects the temperature of the calculated spectrum

        integrated:  `bool`
            if True, plots the integrated/summed spectrum

        saveFile:  `bool`, `str`
            if set, saves the plot to 'saveFile'

        linLog:  `str`
            should be either 'lin' for linear, or 'log' for logarithmic base10

        doLabel: `bool`
            if set to True, labels of the top spectral lines are added

        lw:  `int`, width of the label line in matplotlib units (default=1)

        doTitle:  `bool` if True, then a title is applied to the plot (default=True)


        top:  integer
            specifies to plot only the top strongest lines, default = 10
            if set to None, all lines are plotted

        '''
        fs = 14
        #  fontsize for labels:
        lfs = 12
        plt.figure()
#        mask = self.Em > 1.
#
#        if mask.sum() == 0:
#            self.Ylabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1} \AA^{-1}$ ($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'
#        else:
#            self.Ylabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1} \AA^{-1}$'
        #
#        self.Xlabel = 'Wavelength ('+self.Defaults['wavelength'].capitalize() +')'
        #
        if hasattr(self, 'Spectroscopic'):
            title1 = self.Spectroscopic
        elif hasattr(self, 'ClassName'):
            title1 = self.ClassName
        else:
            title1 = self.__class__.__name__

        if wvlRange is None:
            if hasattr(self, 'WvlRange'):
                wvlRange = self.WvlRange
            else:
                wvlRange = [self.Intensity['wvl'].min(),  self.Intensity['wvl'].max()]
                print(' wvlRange should be specified')


        nTempDens = self.NTempDens




        if integrated:
            lineIntensity = self.Intensity['integrated']
            lineWvl = self.Intensity['wvl']
            lineIonS = self.Intensity['ionS']
            wvlIndex = util.between(lineWvl, wvlRange)
            if hasattr(self, 'Continuum'):
                continuum = self.Continuum['intensity'].sum(axis=0)
            else:
                continuum = np.zeros_like(lineIntensity)

            lineIntensity = lineIntensity[wvlIndex]
            lineWvl = lineWvl[wvlIndex]
            lineIonS = lineIonS[wvlIndex]
        elif nTempDens == 1:
            index = 0
            lineIntensity = self.Intensity['intensity'][0]
            if hasattr(self, 'Continuum'):
                continuum = self.Continuum['intensity']
            else:
                continuum = np.zeros_like(lineIntensity)
            lineWvl = self.Intensity['wvl']
            lineIonS = self.Intensity['ionS']
        else:
            if index is None:
                index = nTempDens//2
            lineIntensity = self.Intensity['intensity'][index]
            lineWvl = self.Intensity['wvl']
            lineIonS = self.Intensity['ionS']
            if hasattr(self, 'Continuum'):
                continuum = self.Continuum['intensity'][index]
            else:
                continuum = np.zeros_like(lineIntensity)
        wvlIndex = util.between(lineWvl, wvlRange)

        lineIntensity = lineIntensity[wvlIndex]
        lineWvl = lineWvl[wvlIndex]
        lineIonS = lineIonS[wvlIndex]
        lineIonSpectr = []
        for anion in lineIonS:
            nameDict = util.convertName(anion)
            lineIonSpectr.append(nameDict['spectroscopic'])
        lineIonSpectr = np.asarray(lineIonSpectr)

        self.Error = 0
        if lineWvl.size == 0:
            print('No lines in this wavelength interval')
            self.Error = 1
            self.Message = 'No lines in this wavelength interval'
            return
        elif top is None:
#            top = lineWvl.size
            top = 0
            lineWvl = []
        elif lineWvl.size > top:
            intsrt = np.argsort(lineIntensity)
            lineWvl = lineWvl[intsrt[-top:]]
            lineIonS = lineIonS[intsrt[-top:]]
            lineIonSpectr = lineIonSpectr[intsrt[-top:]]
            lineIntensity = lineIntensity[intsrt[-top:]]
        else:
            top = lineWvl.size


        #
        plt.ion()
        #
        if integrated:
            spectrum = self.Spectrum['integrated']
            if 'wavelength' in sorted(self.Spectrum.keys()):
                plt.plot(self.Spectrum['wavelength'], spectrum)
            elif 'wvl' in sorted(self.Spectrum.keys()):
                plt.plot(self.Spectrum['wvl'], self.Spectrum['integrated'])
            title = title1 + ' integrated spectrum'
        else:
            nTempDens = self.NTempDens
            if nTempDens == 1:
            #
                spectrum = self.Spectrum['intensity']
                if 'wavelength' in sorted(self.Spectrum.keys()):
                    plt.plot(self.Spectrum['wavelength'], spectrum)
                elif 'wvl' in sorted(self.Spectrum.keys()):
                    plt.plot(self.Spectrum['wvl'], spectrum)
                title = title1 + ' Temperature = %10.2e K'%(self.Temperature)
            else:
                if index is None:
                    index = nTempDens/2
                spectrum = self.Spectrum['intensity'][index]
                if 'wavelength' in sorted(self.Spectrum.keys()):
                    plt.plot(self.Spectrum['wavelength'], spectrum)
                elif 'wvl' in sorted(self.Spectrum.keys()):
                    plt.plot(self.Spectrum['wvl'], spectrum)
                title = title1 + ' Temperature = %10.2e K for index = %3i'%(self.Temperature[index],  index)
        plt.xlabel(self.Spectrum['xlabel'],  fontsize=fs)
        plt.ylabel(self.Spectrum['ylabel'],  fontsize=fs)
        if doTitle:
            plt.title(title, fontsize=fs)

        if doLabel:
            useFilter = self.Spectrum['filter']
            useFactor = self.Spectrum['filterWidth']
            for iwvl, awvl in enumerate(lineWvl):
                filterFactor = useFilter(self.Spectrum['wavelength'], awvl, factor=useFactor).max()
                # need to take into account the continuum
                idx = np.argmin(np.abs(lineWvl[iwvl] - self.Wavelength))
                spIntens = filterFactor*lineIntensity[iwvl] + continuum[idx]
                plt.plot([awvl,  awvl], [0.,  1.2*spIntens], 'k',  lw=lw)
                ypos = 1.25*spIntens
                lbl = lineIonSpectr[iwvl] + ' %8.3f'%(awvl)
                plt.text(awvl,  ypos, lbl, va='bottom', ha='center',rotation='vertical',  fontsize=lfs)

        wdx = util.between(self.Spectrum['wavelength'],  wvlRange)
        ymax = spectrum[wdx].max()
#        plt.xlim(wvlRange)
#        plt.ylim([0., ylim[1]])
        if doLabel:
            plt.axis([wvlRange[0],  wvlRange[1],  0., 1.8*ymax])
        else:
            plt.axis([wvlRange[0],  wvlRange[1],  0., 1.1*ymax])

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
