"""
Base classes used in the `ChiantiPy.core.ion` and `ChiantiPy.core.spectrum`
classes. Mostly printing, plotting and saving routines.
"""
import copy
import time
import sys
import numpy as np

import matplotlib.pyplot as plt
import ChiantiPy.tools.util as util
import ChiantiPy.Gui as chGui
import ChiantiPy.tools.data as chdata


class ionTrails(object):
    """
    Base class for `ChiantiPy.core.ion` and `ChiantiPy.core.spectrum`
    """

    def argCheck(self, temperature=None, eDensity=None, pDensity='default', em = None,  verbose=0):
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


        if em is not None:
            em = np.atleast_1d(em)
            self.Em = em
            if em.size == 1:
                self.Em = np.tile(em,self.NTempDens)

            elif em.size != self.NTempDens:
                raise ValueError('the size of em must be either 1 or the size of the larger of temperature or density %5i'%(self.NTempDens))
        else:
            self.Em = np.ones_like(self.Temperature, np.float64)
        if verbose:
            print('Em size:  %5i'%(self.Em.size))

    def intensityList(self, index=-1, wvlRange=None, wvlRanges=None, top=10, relative=0, outFile=0, rightDigits=4 ):
        """
        List the line intensities. Checks to see if there is an existing Intensity attribute. If it exists, then those values are used.
        Otherwise, the `intensity` method is called.

        This method prints an ASCII table with the following columns:

        1. Ion: the CHIANTI style notation for the ion, e.g. 'c_4' for C IV
        2. lvl1:  the lower level of the transition in the CHIANTI .elvlc file
        3. lvl2:  the upper level of the transition in the CHIANTI .elvlc file
        4. lower:  the notation, usually in LS coupling, of the lower fine
           structure level
        5. upper:  the notation, usually in LS coupling, of the upper fine
           structure level
        6. Wvl(A):  the wavelength of the transition in units as specified in
           the chiantirc file.
        7. Intensity
        8. A value:  the Einstein coefficient for spontaneous emission from
           level 'j' to level 'i'
        9. Obs: indicates whether the CHIANTI database considers this an
           observed line or one obtained from theoretical energy levels

        Regarding the intensity column, if 'flux' in the chiantirc file is set
        to 'energy', the intensity is given by,

        .. math::
           I = \Delta E_{ij}n_jA_{ij}\mathrm{Ab}\\frac{1}{N_e}
           \\frac{N(X^{+m})}{N(X)}\mathrm{EM},

        in units of ergs cm\ :sup:`-2` s\ :sup:`-1` sr\ :sup:`-1`. If 'flux' is set to 'photon',

        .. math::
           I = n_jA_{ij}\mathrm{Ab}\\frac{1}{N_e}\\frac{N(X^{+m})}{N(X)}
           \mathrm{EM},

        where,

        - :math:`\Delta E_{ij}` is the transition energy (ergs)
        - :math:`n_j` is  the fractions of ions in level :math:`j`
        - :math:`A_{ij}` is the Einstein coefficient for spontaneous emission
          from level :math:`j` to level :math:`i` (in s\ :sup:`-1`)
        - :math:`\mathrm{Ab}` is the abundance of the specified element
          relative to hydrogen
        - :math:`N_e` is the electron density (in cm\ :sup:`-3`)
        - :math:`N(X^{+m})/N(X)` is the fractional ionization of ion as a
          function of temperature
        - :math:`\mathrm{EM}` is the emission measure integrated along the
          line-of-sight, :math:`\int\mathrm{d}l\,N_eN_H` (cm\ :sup:`-5`) where
          :math:`N_H` is the density of hydrogen (neutral + ionized)
          (cm\ :sup:`-3`)

        Note that if `relative` is set, the line intensity is relative to the
        strongest line and so the output will be unitless.

        Parameters
        -----------
        index :  `int`,optional
            Index the temperature or eDensity array to use.
            -1 (default) sets the specified value to the middle of the array
        wvlRange : `tuple`
            Wavelength range
        wvlRanges : a tuple, list or array that contains at least 2
            2 element tuples, lists or arrays so that multiple
            wavelength ranges can be specified
        top : `int`
            Number of lines to plot, sorted by descending magnitude.
        relative : `int`
            specifies whether to normalize to strongest line
            default (relative = 0) specified that the intensities should be
            their calculated values
        outFile : `str`
            specifies the file that the intensities should be output to
            default(outFile = 0) intensities are output to the terminal
        rightDigits:  `int`
            specifies the format for the wavelengths for the number of digits
            to right of the decimal place
        """

        if not hasattr(self, 'Intensity'):
            try:
                self.intensity()
            #TODO: specify what exception to catch! or omit the try, catch
            except:
                print(' intensities not calculated and emiss() is unable to calculate them')
                print(' perhaps the temperature and/or eDensity are not set')
                return

        temperature = self.Temperature
        eDensity = self.EDensity
        em = self.Em

        ndens = eDensity.size
        ntemp = temperature.size

        intens = copy.deepcopy(self.Intensity)
        if 'errorMessage' in intens.keys():
            print(' errorMessage:  %s'%(intens['errorMessage']))
            return
        intensity = intens['intensity']
        ionS = intens['ionS']
        wvl = intens['wvl']
        lvl1 = intens['lvl1']
        lvl2 = intens['lvl2']
        pretty1 = intens['pretty1']
        pretty2 = intens['pretty2']
        obs = intens['obs']
        avalue = intens['avalue']

        if ndens == 1 and ntemp == 1:
            if index < 0:
                index = 0
            dstr = ' -  Density = %10.2e (cm$^{-3}$)' %(eDensity[index])
            tstr = ' -  T = %10.2e (K)' %(temperature[index])
            intensity = intensity[index]
            print(' temperature =  %10.2e eDensity = %10.2e'%(temperature[index], eDensity[index]))
        elif ndens == 1 and ntemp > 1:
            if index < 0:
                index = ntemp//2
            print('using index = %5i specifying temperature =  %10.2e'%(index, temperature[index]))
            self.Message = 'using index = %5i specifying temperature =  %10.2e'%(index, temperature[index])

            intensity=intensity[index]
        elif ndens > 1 and ntemp == 1:
            if index < 0:
                index = ntemp//2
            print('using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index]))
            self.Message = 'using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index])
            intensity=intensity[index]
        elif ndens > 1 and ntemp > 1:
            if index < 0:
                index = ntemp//2
            print('using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e em = %10.2e'%(index, temperature[index], eDensity[index], em[index]))
            self.Message = 'using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e = %10.2e'%(index, temperature[index], eDensity[index], em[index])
            intensity=intensity[index]

        if wvlRange:
            wvlIndex=util.between(wvl,wvlRange)
        elif wvlRanges:
            wvlIndex = []
            for awvlRange in wvlRanges:
                wvlIndex.extend(util.between(wvl,awvlRange))
        else:
            wvlIndex = range(wvl.size)

        #  get lines in the specified wavelength range
        intensity = intensity[wvlIndex]
        ionS = ionS[wvlIndex]
        wvl = wvl[wvlIndex]
        lvl1 = lvl1[wvlIndex]
        lvl2 = lvl2[wvlIndex]
        avalue = avalue[wvlIndex]
        pretty1 = pretty1[wvlIndex]
        pretty2 = pretty2[wvlIndex]
        obs = obs[wvlIndex]

        self.Error = 0
        if wvl.size == 0:
            print('No lines in this wavelength interval')
            self.Error = 1
            self.Message = 'No lines in this wavelength interval'
            return

        elif top == 0:
            top = wvl.size
        elif top > wvl.size:
            top = wvl.size

        # sort by intensity
        isrt = np.argsort(intensity)
        ionS = ionS[isrt[-top:]]
        wvl = wvl[isrt[-top:]]
        lvl1 = lvl1[isrt[-top:]]
        lvl2 = lvl2[isrt[-top:]]
        obs = obs[isrt[-top:]]
        intensity = intensity[isrt[-top:]]
        avalue = avalue[isrt[-top:]]
        pretty1 = pretty1[isrt[-top:]]
        pretty2 = pretty2[isrt[-top:]]

        # must follow setting top
        if relative:
            intensity = intensity/intensity[:top].max()

        idx = np.argsort(wvl)
        fmt1 = '%5s %5s %5s %25s - %-25s %12s %12s %12s %3s'
        fmt = '%5s %5i %5i %25s - %-25s %12.' + str(rightDigits) + \
            'f %12.2e %12.2e %1s'
        print('   ')
        print(' ------------------------------------------')
        print('   ')
        print(fmt1%('Ion','lvl1','lvl2','lower','upper','Wvl(A)','Intensity','A value','Obs'))
        for kdx in idx:
            print(fmt%(ionS[kdx], lvl1[kdx], lvl2[kdx], pretty1[kdx], pretty2[kdx], wvl[kdx], intensity[kdx], avalue[kdx], obs[kdx]))
        print('   ')
        print(' ------------------------------------------')
        print('   ')
        #
        self.Intensity['wvlTop'] = wvl[idx]
        self.Intensity['intensityTop'] = intensity[idx]
        if outFile:
            fmt1a = '%5s %5s %5s %25s - %-25s %12s %12s %12s %3s \n'
            fmt = '%5s %5i %5i %25s - %-25s %12.4' + str(rightDigits) + \
            'f %12.2e %12.2e %1s \n'
            outpt = open(outFile, 'w')
            outpt.write(fmt1a%('Ion','lvl1','lvl2','lower','upper','Wvl(A)','Intensity','A value','Obs'))
            for kdx in idx:
                outpt.write(fmt%(ionS[kdx], lvl1[kdx], lvl2[kdx], pretty1[kdx], pretty2[kdx], wvl[kdx], intensity[kdx], avalue[kdx], obs[kdx]))
            outpt.close()


    def intensityPlot(self, index=-1, wvlRange=None, top=10, linLog='lin', relative=False, verbose=False, plotFile=0, em=0 ):
        """
        Plot the line intensities. Uses `Intensity` if it already exists. If
        not, call the `intensity` method.

        Parameters
        -----------
        index:  integer
            specified which value of the temperature array or eDensity array to use.
            default (index=-1) sets the specified value to the middle of the array
        wvlRange:  2 element tuple, list or array determines the wavelength range
        top:  integer
            specifies to plot only the top strongest lines, default = 10
        linLog:  str
            default('lin') produces a plot where the intensity scale is linear
            if set to 'log', produces a plot where the intensity scale is logarithmic
        normalize: = 1 specifies whether to normalize to strongest line
            default (relative = 0) specified that the intensities should be their calculated
            values
        plotFile:
            default=0, the plot is not saved to a file
            othewise, the plot is saved to the 'plotFile'
        em:  emission measure
            if an Intensity attribute needs be created, then the emission measure is applied
        """
        if hasattr(self, 'Spectroscopic'):
            title = self.Spectroscopic
        else:
            title = ''

        if not hasattr(self, 'Intensity'):
            try:
                self.intensity()
            except:
                print(' emissivities not calculated and emiss() is unable to calculate them')
                print(' perhaps the temperature and/or eDensity are not set')
                return

        wvl = self.Intensity['wvl']
        temperature = self.Temperature
        eDensity = self.EDensity
        ndens = eDensity.size
        ntemp = temperature.size

        if ndens == 1 and ntemp == 1:
            if index < 0:
                index = 0
            intensity = self.Intensity['intensity'][index]
            dstr = ' -  Density = %10.2e (cm$^{-3}$)' %(eDensity[index])
            tstr = ' -  T = %10.2e (K)' %(temperature[index])
        elif ndens == 1 and ntemp > 1:
            if index < 0:
                index = ntemp//2
                print('using index = %5i specifying temperature =  %10.2e'%(index, temperature[index]))
                self.Message = 'using index = %5i specifying temperature =  %10.2e'%(index, temperature[index])

            intensity=self.Intensity['intensity'][index]
            dstr=' -  Density = %10.2e (cm$^{-3}$)' % eDensity
            tstr=' -  T = %10.2e (K)' % temperature[index]
        elif ndens > 1 and ntemp == 1:
            if index < 0:
                index = ntemp//2
                print('using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index]))
                self.Message = 'using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index])

            intensity = self.Intensity['intensity'][index]
            dstr=' -  Density = %10.2e (cm$^{-3}$)' % eDensity[index]
            tstr=' -  T = %10.2e (K)' % temperature
        elif ndens > 1 and ntemp > 1:
            if index < 0:
                index = ntemp//2
                print('using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index]))
                self.Message = 'using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index])

            intensity = self.Intensity['intensity'][index]
            dstr=' -  Density = %10.2e (cm$^{-3}$)' % eDensity[index]
            tstr=' -  T = %10.2e (K)' % temperature[index]

        if type(wvlRange) != type(None):
            wvlIndex = util.between(wvl, wvlRange)
        else:
            wvlIndex = range(wvl.size)
        intensity = intensity[wvlIndex]
        wvl = wvl[wvlIndex]

        self.Error = 0
        if wvl.size == 0:
            print('No lines in this wavelength interval')
            self.Error = 1
            self.Message = 'No lines in this wavelength interval'
            return
        elif top == 0:
            top = wvl.size
        elif wvl.size > top:
            intsrt = np.argsort(intensity)
            wvl = wvl[intsrt[-top:]]
            intensity = intensity[intsrt[-top:]]
        else:
            top = wvl.size

        # must follow setting top
        plt.figure()
        ylabel = 'intensity'
        if relative:
            intensity = intensity/intensity[-1]
            ylabel += ' (Relative)'
        xlbl = 'Wavelength ('+self.Defaults['wavelength'] +')'
        ymin = 10.**(np.log10(intensity[0].min()).round(0)-0.5 )
        plt.ion()
        for idx in range(top):
            xx=[wvl[idx], wvl[idx]]
            if linLog == 'lin':
                yy=[0., intensity[idx]]
                plt.plot(xx, yy)
            else:
                yy=[ymin/10., intensity[idx]]
                plt.semilogy(xx, yy)
        plt.xlabel(xlbl)
        plt.ylabel(ylabel)
        plt.title(title+tstr+dstr)
        plt.tight_layout()
        if wvlRange:
            plt.axis([wvlRange[0], wvlRange[1], ymin, intensity.max()])
        if plotFile:
            plt.savefig(plotFile)

    def intensityRatio(self, wvlRange=None, wvlRanges=None, top=10):
        """
        Plot the intensity ratio of 2 lines or sums of lines.
        Shown as a function of density and/or temperature.
        For a single wavelength range, set wvlRange = [wMin, wMax]
        For multiple wavelength ranges, set wvlRanges = [[wMin1,wMax1],[wMin2,wMax2], ...]
        A plot of relative emissivities is shown and then a dialog appears for the user to
        choose a set of lines.

        Parameters
        -----------
        wvlRange : array-like
            Wavelength range, i.e. min and max
        wvlRanges: a tuple, list or array that contains at least 2
            2 element tuples, lists or arrays so that multiple
            wavelength ranges can be specified
        top : `int`
            specifies to plot only the top strongest lines, default = 10
        """

        if not hasattr(self, 'Intensity'):
            try:
                self.intensity()
            except:
                print(' intensities not calculated and emiss() is unable to calculate them')
                print(' perhaps the temperature and/or eDensity are not set')
                return
        fontsize=14
        eDensity = self.EDensity
        temperature = self.Temperature
        ntemp = temperature.size
        if ntemp > 0:
            if temperature[0] == temperature[-1]:
                ntemp = 1
        ndens = eDensity.size
        if ndens > 0:
            if eDensity[0] == eDensity[-1]:
                ndens = 1
        print(' ndens = %5i ntemp = %5i'%(ndens, ntemp))

        ionS = self.Intensity['ionS']
        #  see if we are dealing with more than a single ion
        ionSet = set(ionS)
        ionNum = len(ionSet)
        wvl = self.Intensity["wvl"]
        # find which lines are in the wavelength range if it is set
        if wvlRange:
            igvl=util.between(wvl,wvlRange)
            if len(igvl) == 0:
                print('no lines in wavelength range %12.2f - %12.2f'%(wvlRange[0], wvlRange[1]))
                return
        elif wvlRanges:
            igvl = []
            for awvlRange in wvlRanges:
                igvl.extend(util.between(wvl,awvlRange))
            if len(igvl) == 0:
                print('no lines in wavelength ranges specified ')
                return
        else:
            igvl=range(len(wvl))
        nlines=len(igvl)
        igvl=np.take(igvl,wvl[igvl].argsort())
        if top > nlines:
            top=nlines
        intensity = self.Intensity['intensity']
        maxIntens = np.zeros(nlines,np.float64)
        for iline in range(nlines):
            maxIntens[iline] = intensity[:, igvl[iline]].max()
        for iline in range(nlines):
            if maxIntens[iline]==maxIntens.max():
                maxAll=intensity[:, igvl[iline]]
        igvlsort=np.take(igvl,np.argsort(maxIntens))
        topLines=igvlsort[-top:]
        maxWvl='%5.3f' % wvl[topLines[-1]]
        topLines=topLines[wvl[topLines].argsort()]
        print(' maxWvl = %s'%(maxWvl))

        # need to make sure there are no negative values before plotting
        good = intensity > 0.
        intensMin = intensity[good].min()
        bad = intensity <= 0.
        intensity[bad] = intensMin

        ylabel='Intensity relative to '+maxWvl
        if ionNum == 1:
            title=self.Spectroscopic
        else:
            title = ''

        if ndens==1 and ntemp==1:
            print(' only a single temperature and eDensity')
            return
        elif ndens == 1:
            xlbl='Temperature (K)'
            xvalues=self.Temperature
            outTemperature=self.Temperature
            outDensity = self.EDensity
            desc_str=' Density = %10.2e (cm$^{-3}$)' % self.EDensity[0]
        elif ntemp == 1:
            xvalues=self.EDensity
            outTemperature = self.Temperature
            outDensity=self.EDensity
            xlbl='Electron Density (cm$^{-3}$)'
            desc_str=' Temp = %10.2e (K)' % self.Temperature[0]
        else:
            outTemperature=self.Temperature
            outDensity=self.EDensity
            xlbl='Temperature (K)'
            xvalues=self.Temperature
            desc_str=' Variable Density'
        # put all actual plotting here
        plt.ion()
        #  maxAll is an array
        ymax = np.max(intensity[:, topLines[0]]/maxAll)
        ymin = ymax
        plt.figure()
        ax = plt.subplot(111)
        nxvalues=len(xvalues)
        # reversing is necessary - otherwise, get a ymin=ymax and a matplotlib error
        for iline in range(top-1, -1, -1):
            tline=topLines[iline]
            plt.loglog(xvalues,intensity[:, tline]/maxAll)
            if np.min(intensity[:, tline]/maxAll) < ymin:
                ymin = np.min(intensity[:, tline]/maxAll)
            if np.max(intensity[:, tline]/maxAll) > ymax:
                ymax = np.max(intensity[:, tline]/maxAll)
            skip=2
            start=divmod(iline,nxvalues)[1]
            for ixvalue in range(start,nxvalues,nxvalues//skip):
                if ionNum == 1:
                    text = '%10.4f'%(wvl[tline])
                else:
                    text = '%s %10.4f'%(ionS[tline], wvl[tline])
                plt.text(xvalues[ixvalue], intensity[ixvalue, tline]/maxAll[ixvalue], text)
        if ndens == 1:
            print('%12.2e  %12.2e '%(xvalues.min(),xvalues.max()))
            plt.xlim(xvalues.min(),xvalues.max())
            plt.xlabel(xlbl,fontsize=fontsize)
            plt.ylabel(ylabel,fontsize=fontsize)
        elif ntemp == 1:
            plt.text(0.07, 0.5,title, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            ax2 = plt.twiny()
            xlblDen=r'Electron Density (cm$^{-3}$)'
            plt.xlabel(xlblDen, fontsize=fontsize)
            plt.loglog(eDensity,intensity[:, topLines[top-1]]/maxAll, visible=False)
            ax2.xaxis.tick_top()
            plt.ylim(ymin/1.2, 1.2*ymax)
        else:
            plt.ylim(ymin/1.2, 1.2*ymax)
            plt.title(title+desc_str,fontsize=fontsize)
        plt.tight_layout()
        plt.draw()
        #  need time to let matplotlib finish plotting
        time.sleep(0.5)

        # get line selection
        selectTags = []
        for itop in topLines:
            if ionNum == 1:
                selectTags.append(str(wvl[itop]))
            else:
                selectTags.append(ionS[itop]+ ' '+ str(wvl[itop]))
        numden = chGui.gui.choice2Dialog(selectTags)

        # num_idx and den_idx are tuples
        num_idx=numden.numIndex
        if len(num_idx) == 0:
            print(' no numerator lines were selected')
            return
        den_idx=numden.denIndex
        if len(den_idx) == 0:
            print(' no denominator lines were selected')
            return
        numIntens=np.zeros(len(xvalues),np.float64)
        for aline in num_idx:
            numIntens += intensity[:, topLines[aline]]

        denIntens = np.zeros(len(xvalues),np.float64)
        for aline in den_idx:
            denIntens += intensity[:, topLines[aline]]

        # plot the desired ratio
        #  maxAll is an array
        plt.figure()
        ax = plt.subplot(111)
        plt.loglog(xvalues,numIntens/denIntens)
        plt.xlim(xvalues.min(),xvalues.max())
        plt.xlabel(xlbl,fontsize=fontsize)
        plt.ylabel('Ratio ('+self.Defaults['flux']+')',fontsize=fontsize)
        if ionNum == 1:
            desc = ionS[0]
        else:
            desc = ''
        for aline in num_idx:
            if ionNum == 1:
                desc += ' ' + str(wvl[topLines[aline]])
            else:
                desc += ' ' + ionS[topLines[aline]] + ' ' + str(wvl[topLines[aline]])
        desc += ' / '
        for aline in den_idx:
            if ionNum == 1:
                desc += ' ' + str(wvl[topLines[aline]])
            else:
                desc += ' ' + ionS[topLines[aline]] + ' ' + str(wvl[topLines[aline]])
        if ndens == ntemp and ntemp > 1:
            plt.text(0.07, 0.5,desc, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            #
            ax2 = plt.twiny()
            xlblDen=r'Electron Density (cm$^{-3}$)'
            plt.xlabel(xlblDen, fontsize=fontsize)
            plt.loglog(eDensity,numIntens/denIntens, visible=False)
            ax2.xaxis.tick_top()
        else:
            plt.title(desc,fontsize=fontsize)
        plt.tight_layout()

        cnt = desc.count(' ')
        intensityRatioFileName = desc.replace(' ', '_', cnt) + '.rat'
        intensityRatioFileName = intensityRatioFileName.lstrip('_').replace('_/_','-')
        self.IntensityRatio={'ratio':numIntens/denIntens,'desc':desc,
                'temperature':outTemperature,'eDensity':outDensity,'filename':intensityRatioFileName, 'numIdx':num_idx, 'denIdx':den_idx}

    def intensityRatioSave(self,outFile=0):
        """
        Save the intensity ratio to a file.

        The intensity ratio as a function to temperature and eDensity is saved to an asciii file. Descriptive information is included at the top of the file.

        Parameters
        ----------
            outFile
                default(0):  the plot of the intensity ratio is not saved
                str/unicode:  the plot is saved to the file names 'outFile'

        """
        if not outFile:
            outFile = self.IntensityRatio['filename']
            print(' saving ratio to filename = %s'%(outFile))
        if hasattr(self, 'IntensityRatio'):
            temperature=self.IntensityRatio['temperature']
            eDensity=self.IntensityRatio['eDensity']
            ratio=self.IntensityRatio['ratio']
            out = open(outFile,'w')
            nvalues=len(ratio)
            #  need to add 7 lines to maintain IDL like files
            out.write(outFile+'\n')    #1
            out.write(self.IntensityRatio['desc']+'\n') #2
            out.write(' created with ChiantiPy using CHIANTI version '+ chdata.ChiantiVersion +'\n')   #3
            out.write(' columns are temperature, eDensity, ratio'+'\n')  #5
            tunit = 'K'
            out.write(' temperature in '+tunit+', electron eDensity in cm^(-3)'+'\n')  #6
            out.write(' ratio given in '+self.Defaults['flux']+'\n')   #4
            out.write(' '+'\n') #7
            for ivalue in range(nvalues):
                s='%12.3e %12.3e  %12.3e \n' % (temperature[ivalue],eDensity[ivalue],ratio[ivalue])
                out.write(s)
            out.close()
        else:
            print(' in .intensityRatioSave(), no IntensityRatio is found')

