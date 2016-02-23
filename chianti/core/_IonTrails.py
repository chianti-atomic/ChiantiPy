'''
_IonTrails
'''
import copy
import time
import numpy as np
import pylab as pl
import chianti.util as util
import chianti.io as io
import chianti.Gui as chGui
import chianti.data as chdata
#
class _ionTrails():
    '''
    a collection of methods for use in ion and spectrum calculations
    '''
    def __init__(self):
        pass
        return
        #
        # -------------------------------------------------------------------------------------
        #
#    def gofnt(self,wvlRange=0,top=10, verbose=0):
#        """
#        Calculate the 'so-called' G(T) function.
#
#        Given as a function of both temperature and eDensity.
#
#        Only the top( set by 'top') brightest lines are plotted.
#        the G(T) function is returned in a dictionary self.Gofnt
#        """
#        #
#        #self.emiss={"wvl":wvl,"emiss":em,"units":units,"plotLabels":plotLabels}
#        #
#        #
#        #
#        if not hasattr(self, 'Intensity'):
#            try:
#                self.intensity()
#            except:
#                print(' intensities not calculated and emiss() is unable to calculate them')
#                print(' perhaps the temperature and/or eDensity are not set')
#                return
#        #
#        # everything in self.Intensity should be a numpy array
#        #
#         #
#            #
#        #
#        #
##        if hasattr(self, 'Abundance'):
##            ab=self.Abundance
##        else:
##            self.Abundance = io.abundanceRead()
##            ab=self.Abundance
##        if not hasattr(self, 'Abundance'):
##            self.Abundance = io.abundanceRead()
#        #
#        fontsize=12
#        #
##        emiss=em["emiss"]
##        wvl=em["wvl"]
##        pretty1 = em['pretty1']
##        pretty2 = em['pretty2']
##        lvl1 = em['lvl1']
##        lvl2 = em['lvl2']
#        #
#        intens = copy.deepcopy(self.Intensity)
#        intensity = intens['intensity']
#        ionS = intens['ionS']        
#        ionSet = set(ionS)
#        ionNum = len(ionSet)
#        wvl = intens['wvl']
#        lvl1 = intens['lvl1']
#        lvl2 = intens['lvl2']
#        pretty1 = intens['pretty1']
#        pretty2 = intens['pretty2']
#        obs = intens['obs']
#        avalue = intens['avalue']
#        #
#        temperature=self.Temperature
#        eDensity=self.EDensity
##        plotLabels=em["plotLabels"]
##        xLabel=plotLabels["xLabel"]
##        yLabel=plotLabels["yLabel"]
#        #
#        # find which lines are in the wavelength range if it is set
#        #
#        #
#        if type(wvlRange) != type(1):
#            igvl=util.between(wvl,wvlRange)
#        else:
#            igvl=range(len(wvl))
#        nlines=len(igvl)
#        if nlines ==0:
#            print(' no lines in selected interval')
#            return
#        # find the top most intense lines
#        #
#        if top > nlines:
#            top=nlines
#        maxEmiss=np.zeros(nlines,'Float64')
#        for iline in range(nlines):
#            maxEmiss[iline] = intensity[:, igvl[iline]].max()
#        for iline in range(nlines):
#            if maxEmiss[iline] >= maxEmiss.max():
#                maxAll = intensity[:, igvl[iline]]
##                maxIndex = igvl[iline]
##        print ' maxIndex, maxAll = ', maxIndex,  maxAll
##        line=range(nlines)
#        igvlsort=np.take(igvl,np.argsort(maxEmiss))
#        topLines=igvlsort[-top:]
#        maxWvl='%5.3f' % wvl[topLines[-1]]
##        maxline=topLines[-1]
#        #
#        # need to make sure there are no negative values before plotting
#        good = np.where(intensity > 0.)
#        emissMin = intensity[good].min()
#        bad = np.where(intensity <= 0.)
#        emiss[bad] = emissMin
#        #
#        topLines=topLines[wvl[topLines].argsort()]
#        #
#        #
#        ntemp=self.Temperature.size
#        #
#        ndens=self.EDensity.size
#        #
#        ylabel = 'Emissivity relative to '+maxWvl
#        title = self.Spectroscopic
#        #
#        #
#        if ndens==1 and ntemp==1:
#            print(' only a single temperature and eDensity')
#            return
#        elif ndens == 1:
#            xlabel='Temperature (K)'
#            ngofnt = temperature.size
#            xvalues=temperature
#            outTemperature=temperature
#            outDensity=np.zeros(ntemp,'Float64')
#            outDensity.fill(eDensity)
#            desc_str=' at Density = %10.2e' % eDensity
#        elif ntemp == 1:
#            xvalues=eDensity
#            ngofnt = eDensity.size
#            outTemperature=np.zeros(ndens,'Float64')
#            outTemperature.fill(temperature)
#            outDensity=eDensity
#            xlabel=r'$\rm{Electron Density (cm}^{-3}\rm{)}$'
#            desc_str=' at Temperature = %10.2e' % temperature
#        else:
#            outTemperature=temperature
#            outDensity=eDensity
#            xlabel='Temperature (K)'
#            xvalues=temperature
#            ngofnt = temperature.size
#            desc_str=' for variable Density'
#            #
#        #
#        # put all actual plotting here
#        #
#        pl.ion()
##        if chInteractive:
##            pl.ion()
##        else:
##            pl.ioff()
#        #
#        pl.figure()
#        ax = pl.subplot(111)
#        nxvalues=len(xvalues)
#        #  maxAll is an array
##        print ' emiss = ', np.max(emiss[top-1]), np.max(emiss[0])
##        print ' maxAll = ', maxAll
##        ymax = np.max(1.2*emiss[top-1]/maxAll)
#        ymax = 1.2
##        print ' ymax = ', ymax
#        ymin = ymax  #  np.min(emiss[0]/maxAll)  # was originally  = ymax
#        for iline in range(top):
#            tline=topLines[iline]
#            pl.loglog(xvalues,intensity[tline]/maxAll)
#            if np.min(intensity[:, tline]/maxAll) < ymin:
#                ymin = np.min(intensity[:, tline]/maxAll)
#            skip=2
#            start=divmod(iline,nxvalues)[1]
#            for ixvalue in range(start,nxvalues,nxvalues/skip):
#                pl.text(xvalues[ixvalue],intensity[:, tline,ixvalue]/maxAll[ixvalue],str(wvl[tline]))
#        pl.xlim(xvalues.min(),xvalues.max())
#        pl.ylim(ymin, ymax)
##       yl=pl.ylim()
##       pl.ylim(yl[0],1.2)
#        pl.xlabel(xlabel,fontsize=fontsize)
#        pl.ylabel(ylabel,fontsize=fontsize)
#        if ndens == ntemp and ntemp > 1:
#            pl.text(0.07, 0.5,title, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
#            #
#            ax2 = pl.twiny()
#            xlabelDen=r'Electron Density (cm$^{-3}$)'
#            pl.xlabel(xlabelDen, fontsize=fontsize)
#            pl.loglog(eDensity,intensity[:, topLines[top-1]]/maxAll, visible=False)
#            ax2.xaxis.tick_top()
#        else:
#            pl.ylim(ymin, ymax)
#            pl.title(title+desc_str,fontsize=fontsize)
#        pl.draw()
#        #
#        time.sleep(0.5)
#        #
##        print ' topLInes = ', wvl[topLines]
##        wvlChoices = []
##        for iline in range(top):
##            tline = topLines[iline]
##            wvlChoices.append('%12.4f %4i %4i %s - %s'%(wvl[tline], lvl1[tline], lvl2[tline], pretty1[tline], pretty2[tline]))
##        gline = chGui.gui.selectorDialog(wvlChoices,label='Select line(s)')
#        #
#        selectTags = []
#        for itop in topLines:
#            if ionNum == 1:
#                selectTags.append(str(wvl[itop]))
#            else:
#                selectTags.append(ionS[itop]+ ' '+ str(wvl[itop]))
#        gline = chGui.gui.selectorDialog(selectTags,label='Select line(s)')
#
#        gline_idx=gline.selectedIndex
#        #
#        #
##        gAbund=self.Abundance
#        #
##        try:
##            thisIoneq=self.IoneqOne
##        except:
##            self.ioneqOne()
##            thisIoneq=self.IoneqOne
##        if verbose:
##            print(' abundance = %10.2e '%(gAbund))
##            print(' index  temperature  ion fraction')
##            for it,  anioneq in enumerate(thisIoneq):
##                print (' %5i %10.2e %10.2e '%(it, outTemperature[it], anioneq))
#        #        gioneq=np.where(thisIoneq > 0.)
#        #        y2=interpolate.splrep(np.log(self.IoneqAll['ioneqTemperature'][gioneq]),np.log(thisIoneq[gioneq]),s=0)  #allow smoothing,s=0)
#        #        gIoneq=interpolate.splev(np.log(temperature),y2)   #,der=0)
##        gIoneq=self.IoneqOne/eDensity
#        #
#        #
#        #
#        # plot the desired ratio
#        pl.figure()
#        g_line = topLines[gline_idx]#  [0]
#        ##        print ' g_line = ',g_line
#        #
#        gofnt=np.zeros(ngofnt,'float64')
#        for aline in g_line:
#            gofnt += intensity[:, aline].squeeze()
#            gofnt += intensity[:, aline].squeeze()
#        self.Gofnt={'temperature':outTemperature,'eDensity':outDensity,'gofnt':gofnt, 'index':g_line, 'wvl':wvl[g_line]}
#        #
#        pl.loglog(xvalues,gofnt)
#        pl.xlim(xvalues.min(),xvalues.max())
#        pl.xlabel(xlabel,fontsize=fontsize)
#        pl.ylabel('Gofnt',fontsize=fontsize)
#        if ionNum == 1:
#            newTitle = '%9s'%(self.Spectroscopic) + '%12.3f %4i %4i %s - %s'%(wvl[g_line[0]], lvl1[g_line[0]], lvl2[g_line[0]], pretty1[g_line[0]], pretty2[g_line[0]])
#            if len(g_line) > 1:
#                newTitle +='\n'
#            for igl in g_line[1:]:
#                newTitle += ' ' + '%12.3f %4i %4i %s - %s'%(wvl[igl], lvl1[igl], lvl2[igl], pretty1[igl], pretty2[igl])
#                if igl != g_line[-1]:
#                    newTitle +='\n'
#    #        pl.annotate(newTitle, xytext=(0.3, 0.3), textcoords='figure_fraction')
#        else:
#            newTitle = 'newTitle'
#        pl.annotate(newTitle, xy=(-10, 10),
#                xycoords='axes points',
#                horizontalalignment='right', verticalalignment='bottom')  #,fontsize=20)
#        if ndens == ntemp and ntemp > 1:
##            newTitle +' '+str(wvl[g_line])+' '+desc_str
#            pl.text(0.07, 0.5,newTitle, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
#            #
#            ax2 = pl.twiny()
##            xlabel=r'Electron Density (cm$^{-3}$)'
#            pl.xlabel(xlabelDen, fontsize=fontsize)
#            pl.loglog(eDensity,gofnt, visible=False)
#            ax2.xaxis.tick_top()
#        else:
#            pl.title(newTitle, fontsize=fontsize)
#        #pl.ioff()
#        #pl.show()
##        return
        #
        # ---------------------------------------------------------------------------
        #
    def intensityList(self, index=-1,  wvlRange=None, wvlRanges=None,   top=10, relative=0, outFile=0 ):
        '''
        List the line intensities

        wvlRange, a 2 element tuple, list or array determines the wavelength range

        Top specifies to plot only the top strongest lines, default = 10

        normalize = 1 specifies whether to normalize to strongest line, default = 0
        rewrite of emissList
        '''
        #
        #
        if not hasattr(self, 'Intensity'):
            try:
                self.intensity()
            except:
                print(' intensities not calculated and emiss() is unable to calculate them')
                print(' perhaps the temperature and/or eDensity are not set')
                return
        #
        # everything in self.Intensity should be a numpy array
        #
         #
            #
        temperature = self.Temperature
        eDensity = self.EDensity
        em = self.Em
        #
        ndens = eDensity.size
        ntemp = temperature.size
        #
        #
        intens = copy.deepcopy(self.Intensity)
        intensity = intens['intensity']
        ionS = intens['ionS']
        wvl = intens['wvl']
        lvl1 = intens['lvl1']
        lvl2 = intens['lvl2']
        pretty1 = intens['pretty1']
        pretty2 = intens['pretty2']
        obs = intens['obs']
        avalue = intens['avalue']
        #
        if ndens == 1 and ntemp == 1:
            dstr = ' -  Density = %10.2e (cm$^{-3}$)' %(eDensity)
            tstr = ' -  T = %10.2e (K)' %(temperature)
        elif ndens == 1 and ntemp > 1:
            if index < 0:
                index = ntemp/2
            print('using index = %5i specifying temperature =  %10.2e'%(index, temperature[index]))
            self.Message = 'using index = %5i specifying temperature =  %10.2e'%(index, temperature[index])

            intensity=intensity[index]
        elif ndens > 1 and ntemp == 1:
            if index < 0:
                index = ntemp/2
            print('using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index]))
            self.Message = 'using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index])
            intensity=intensity[index]
        elif ndens > 1 and ntemp > 1:
            if index < 0:
                index = ntemp/2
            print('using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e em = %10.2e'%(index, temperature[index], eDensity[index], em[index]))
            self.Message = 'using index = %5i specifying temperature = %10.2e, eDensity =  %10.2eem = %10.2e'%(index, temperature[index], eDensity[index], em[index])
            intensity=intensity[index]
        #
#        print('shpae of intensity 0 = %5i %5i'%(intensity.shape))
        #
        if wvlRange:
            wvlIndex=util.between(wvl,wvlRange)
        elif wvlRanges:
            wvlIndex = []
            for awvlRange in wvlRanges:
                wvlIndex.extend(util.between(wvl,awvlRange))
        else:
            wvlIndex = range(wvl.size)
        #
        #  get lines in the specified wavelength range
        #
#        print('shpae of intensity 1 = %5i %5i'%(intensity.shape))
        #
        intensity = intensity[wvlIndex]
        ionS = ionS[wvlIndex]
        wvl = wvl[wvlIndex]
        lvl1 = lvl1[wvlIndex]
        lvl2 = lvl2[wvlIndex]
        avalue = avalue[wvlIndex]
        pretty1 = pretty1[wvlIndex]
        pretty2 = pretty2[wvlIndex]
        obs = obs[wvlIndex]
        #
        self.Error = 0
        if wvl.size == 0:
            print('No lines in this wavelength interval')
            self.Error = 1
            self.Message = 'No lines in this wavelength interval'
            return
        #
        elif top == 0:
            top = wvl.size
        elif top > wvl.size:
            top = wvl.size
        #
        # sort by intensity
        #
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
        #
#        print('shpae of intensity 2= %5i %5i'%(intensity.shape))
        #
    # must follow setting top
        #
        if relative:
            intensity = intensity/intensity[:top].max()
        #
        #
        idx = np.argsort(wvl)
        fmt1 = '%5s %5s %5s %25s - %-25s %12s %12s %12s %3s'
        fmt = '%5s %5i %5i %25s - %-25s %12.4f %12.3e %12.2e %1s'
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
            fmt = '%5s %5i %5i %25s - %-25s %12.4f %12.3e %12.2e %1s \n'
            outpt = open(outFile, 'w')
            outpt.write(fmt1a%('Ion','lvl1','lvl2','lower','upper','Wvl(A)','Intensity','A value','Obs'))
            for kdx in idx:
                outpt.write(fmt%(ionS[kdx], lvl1[kdx], lvl2[kdx], pretty1[kdx], pretty2[kdx], wvl[kdx], intensity[kdx], avalue[kdx], obs[kdx]))
            outpt.close()
        #
        # ---------------------------------------------------------------------------
        #
    def intensityPlot(self, index=-1,  wvlRange=None,  top=10, linLog='lin', relative=0,  verbose=0, plotFile = 0, saveFile=0, em=0 ):
        '''Plot the line intensities.

        wvlRange, a 2 element tuple, list or array determines the wavelength range

        Top specifies to plot only the top strongest lines, default = 10

        linLog specifies a linear or log plot, want either lin or log, default = lin

        normalize = 1 specifies whether to normalize to strongest line, default = 0'''
        #
        title=self.Spectroscopic
        #
        if hasattr(self, 'Intensity'):
            intens = self.Intensity['intensity']
        else:
            try:
                self.intensity(em=em)
                intens = self.Intensity['intensity']
            except:
                print(' emissivities not calculated and emiss() is unable to calculate them')
                print(' perhaps the temperature and/or eDensity are not set')
                return
        wvl = self.Intensity['wvl']
        temperature = self.Temperature
        eDensity = self.EDensity
        #
        ndens = eDensity.size
        ntemp = temperature.size
        #
        if ndens == 1 and ntemp == 1:
            dstr = ' -  Density = %10.2e (cm$^{-3}$)' %(eDensity)
            tstr = ' -  T = %10.2e (K)' %(temperature)
        elif ndens == 1 and ntemp > 1:
            if index < 0:
                index = ntemp/2
                print('using index = %5i specifying temperature =  %10.2e'%(index, temperature[index]))
                self.Message = 'using index = %5i specifying temperature =  %10.2e'%(index, temperature[index])
#            if chInteractive:
#                print 'using index = %5i specifying temperature =  %10.2e'%(index, temperature[index])
#            else:
#                self.Message = 'using index = %5i specifying temperature =  %10.2e'%(index, temperature[index])
            intensity=self.Intensity['intensity'][index]
            dstr=' -  Density = %10.2e (cm$^{-3}$)' % eDensity
            tstr=' -  T = %10.2e (K)' % temperature[index]
        elif ndens > 1 and ntemp == 1:
            if index < 0:
                index = ntemp/2
                print('using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index]))
                self.Message = 'using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index])
#            if chInteractive:
#                print 'using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index])
#            else:
#                self.Message = 'using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index])
            intensity = self.Intensity['intensity'][index]
            dstr=' -  Density = %10.2e (cm$^{-3}$)' % eDensity[index]
            tstr=' -  T = %10.2e (K)' % temperature
        elif ndens > 1 and ntemp > 1:
            if index < 0:
                index = ntemp/2
                print('using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index]))
                self.Message = 'using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index])
#             if chInteractive:
#                print 'using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index])
#            else:
#                self.Message = 'using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index])
            intensity = self.Intensity['intensity'][index]
            dstr=' -  Density = %10.2e (cm$^{-3}$)' % eDensity[index]
            tstr=' -  T = %10.2e (K)' % temperature[index]
        if type(wvlRange) != type(None):
            wvlIndex = util.between(wvl, wvlRange)
        else:
            wvlIndex = range(wvl.size)
        intensity = intensity[wvlIndex]
        wvl = wvl[wvlIndex]
        #
        self.Error = 0
        if wvl.size == 0:
            print('No lines in this wavelength interval')
            self.Error = 1
            self.Message = 'No lines in this wavelength interval'
#            if chInteractive:
#                print 'No lines in this wavelength interval'
#            else:
#                self.Error = 1
#                self.Message = 'No lines in this wavelength interval'
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
        #
        pl.figure()
        ylabel = 'intensity'
        if relative:
            intensity = intensity/intensity[-1]
            ylabel += ' (Relative)'
        #
        xlabel = 'Wavelength ('+self.Defaults['wavelength'] +')'
        #
        ymin = 10.**(np.log10(intensity[0].min()).round(0)-0.5 )
        #
        pl.ion()
#        if chInteractive:
#            pl.ion()
#        else:
#            pl.ioff()
        #
        for idx in range(top):
            xx=[wvl[idx], wvl[idx]]
            if linLog == 'lin':
                yy=[0., intensity[idx]]
                pl.plot(xx, yy)
            else:
                yy=[ymin/10., intensity[idx]]
                pl.semilogy(xx, yy)
        pl.xlabel(xlabel)
        pl.ylabel(ylabel)
        pl.title(title+tstr+dstr)
        if wvlRange:
            pl.axis([wvlRange[0], wvlRange[1], ymin, intensity.max()])
        if plotFile:
            pl.savefig(plotFile)
        #
#        idx = np.argsort(wvl)
#        self.Intens['wvlTop'] = wvl[idx]
#        self.Intensity['intensTop'] = intensity[idx]
        #
        # -------------------------------------------------------------------------------------
        #
    def intensityRatio(self,wvlRange=None, wvlRanges=None,top=10):
        """
        Plot the ratio of 2 lines or sums of lines.
        Shown as a function of density and/or temperature.
        For a single wavelength range, set wvlRange = [wMin, wMax]
        For multiple wavelength ranges, set wvlRanges = [[wMin1,wMax1],[wMin2,wMax2], ...]
        A plot of relative emissivities is shown and then a dialog appears for the user to
        choose a set of lines.
        """
        #
        #        self.Emiss={"temperature":temperature,"density":density,"wvl":wvl,"emiss":em,
        #        "plotLabels":plotLabels}
        #
        #
        if not hasattr(self, 'Intensity'):
            try:
                self.intensity()
            except:
                print(' intensities not calculated and emiss() is unable to calculate them')
                print(' perhaps the temperature and/or eDensity are not set')
                return
        #
        # everything in self.Intensity should be a numpy array
        #
#        intens = copy.copy(self.Intensity)
#        intensity = intens['intensity']
        #
        #
        fontsize=14
        #
#        temperature = self.Temperature
        eDensity = self.EDensity
        temperature = self.Temperature
#        intensity = intens['intensity']
        #
#        temp=np.asarray(temperature,'Float32')
        ntemp = temperature.size
        if ntemp > 0:
            if temperature[0] == temperature[-1]:
                ntemp = 1
        #
        ndens = eDensity.size
        if ndens > 0:
            if eDensity[0] == eDensity[-1]:
                ndens = 1
        #
        print(' ndens = %5i ntemp = %5i'%(ndens, ntemp))
        
        ionS = self.Intensity['ionS']
        #  see if we are dealing with more than a single ion
        ionSet = set(ionS)
        ionNum = len(ionSet)
        wvl = self.Intensity["wvl"]
#        plotLabels = intens["plotLabels"]
#        xLabel = plotLabels["xLabel"]
#        yLabel = plotLabels["yLabel"]
        #
        # find which lines are in the wavelength range if it is set
        #
        #
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
        #
        nlines=len(igvl)
        #
#        print ' nlines = ',nlines
#        print ' iglv = ',igvl
        igvl=np.take(igvl,wvl[igvl].argsort())
        # find the top most intense lines
        #
        if top > nlines:
            top=nlines
            #
        intensity = self.Intensity['intensity']
        maxIntens = np.zeros(nlines,'Float64')
        for iline in range(nlines):
            maxIntens[iline] = intensity[:, igvl[iline]].max()
        for iline in range(nlines):
            if maxIntens[iline]==maxIntens.max():
                maxAll=intensity[:, igvl[iline]]
#        line=range(nlines)
        igvlsort=np.take(igvl,np.argsort(maxIntens))
#        print 'igvlsort = ', igvlsort
        topLines=igvlsort[-top:]
#        print ' topLines = ', topLines
        maxWvl='%5.3f' % wvl[topLines[-1]]
#        maxline=topLines[-1]
        #
        topLines=topLines[wvl[topLines].argsort()]
        #
        #
        # need to make sure there are no negative values before plotting
        good = intensity > 0.
        intensMin = intensity[good].min()
        bad = intensity <= 0.
        intensity[bad] = intensMin
        #
        #
#        ntemp=self.Temperature.size
        #
#        ndens=self.EDensity.size
        #
        ylabel='Intensity relative to '+maxWvl
        if ionNum == 1:
            title=self.Spectroscopic
        else:
            title = ''
        #
        #
        if ndens==1 and ntemp==1:
            print(' only a single temperature and eDensity')
            return
        elif ndens == 1:
            xlabel='Temperature (K)'
            xvalues=self.Temperature
            outTemperature=self.Temperature
#            outDensity=np.zeros(ntemp,'Float64')
#            outDensity.fill(self.EDensity[0])
            outDensity = self.EDensity
            desc_str=' Density = %10.2e (cm)$^{-3}$' % self.EDensity[0]
        elif ntemp == 1:
            xvalues=self.EDensity
#            outTemperature=np.zeros(ndens,'Float64')
#            outTemperature.fill(self.Temperature[0])
            outTemperature = self.Temperature
            outDensity=self.EDensity
            xlabel=r'$\rm{Electron Density (cm)^{-3}}$'
            desc_str=' Temp = %10.2e (K)' % self.Temperature[0]
        else:
            outTemperature=self.Temperature
            outDensity=self.EDensity
            xlabel='Temperature (K)'
            xvalues=self.Temperature
            desc_str=' Variable Density'
        #
        # put all actual plotting here
        #
        pl.ion()
#        if chInteractive:
#            pl.ion()
#        else:
#            pl.ioff()
        #
        #  maxAll is an array
        ymax = np.max(intensity[:, topLines[0]]/maxAll)
        ymin = ymax
        pl.figure()
        ax = pl.subplot(111)
        nxvalues=len(xvalues)
        for iline in range(top):
            tline=topLines[iline]
            pl.loglog(xvalues,intensity[:, tline]/maxAll)
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
                pl.text(xvalues[ixvalue], intensity[ixvalue, tline]/maxAll[ixvalue], text)
        pl.xlim(xvalues.min(),xvalues.max())
#        pl.ylim(ymin, ymax)
        pl.xlabel(xlabel,fontsize=fontsize)
        pl.ylabel(ylabel,fontsize=fontsize)
        if ndens == ntemp and ntemp > 1:
            pl.text(0.07, 0.5,title, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            #
            ax2 = pl.twiny()
            xlabelDen=r'Electron Density (cm$^{-3}$)'
            pl.xlabel(xlabelDen, fontsize=fontsize)
            pl.loglog(eDensity,intensity[:, topLines[top-1]]/maxAll, visible=False)
            ax2.xaxis.tick_top()
            pl.ylim(ymin/1.2, 1.2*ymax)
        else:
            pl.ylim(ymin/1.2, 1.2*ymax)
            pl.title(title+desc_str,fontsize=fontsize)
        pl.draw()
        #  need time to let matplotlib finish plotting
        time.sleep(0.5)
        #
        # get line selection
        #
        selectTags = []
        for itop in topLines:
            if ionNum == 1:
                selectTags.append(str(wvl[itop]))
            else:
                selectTags.append(ionS[itop]+ ' '+ str(wvl[itop]))
        #
#        numden = chGui.gui.choice2Dialog(wvl[topLines])
        numden = chGui.gui.choice2Dialog(selectTags)
        #
        # num_idx and den_idx are tuples
        #
        num_idx=numden.numIndex
        if len(num_idx) == 0:
            print(' no numerator lines were selected')
            return
        #
        den_idx=numden.denIndex
        if len(den_idx) == 0:
            print(' no denominator lines were selected')
            return
        #
        numIntens=np.zeros(len(xvalues),'Float64')
        for aline in num_idx:
            numIntens += intensity[:, topLines[aline]]
        #
        denIntens = np.zeros(len(xvalues),'Float64')
        for aline in den_idx:
            denIntens += intensity[:, topLines[aline]]
        #
        # plot the desired ratio
        #  maxAll is an array
        pl.figure()
        ax = pl.subplot(111)
        pl.loglog(xvalues,numIntens/denIntens)
        pl.xlim(xvalues.min(),xvalues.max())
        pl.xlabel(xlabel,fontsize=fontsize)
        pl.ylabel('Ratio ('+self.Defaults['flux']+')',fontsize=fontsize)
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
            pl.text(0.07, 0.5,desc, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            #
            ax2 = pl.twiny()
            xlabelDen=r'Electron Density (cm$^{-3}$)'
            pl.xlabel(xlabelDen, fontsize=fontsize)
            pl.loglog(eDensity,numIntens/denIntens, visible=False)
            ax2.xaxis.tick_top()
        else:
#            pl.ylim(ymin, ymax)
            pl.title(desc,fontsize=fontsize)
#       desc=title+' '+str(wvl[num_line])+' / '+str(wvl[den_line])+' '+desc_str
#        pl.title(desc, fontsize=fontsize)
#       pl.title(title+' '+str(wvl[num_line])+' / '+str(wvl[den_line])+' '+desc_str,fontsize=fontsize)
#        pl.draw()
#        pl.ioff()
#        pl.show()
        #
#        intensityRatioFileName=self.IonStr
#        for aline in num_idx:
#            intensityRatioFileName+= '_%3i'%(wvl[topLines[aline]])
#        intensityRatioFileName+='_2'
#        for aline in den_idx:
#            intensityRatioFileName+= '_%3i'%(wvl[topLines[aline]])
        cnt = desc.count(' ')
        intensityRatioFileName = desc.replace(' ', '_', cnt) + '.rat'
        intensityRatioFileName = intensityRatioFileName.lstrip('_').replace('_/_','-')
        self.IntensityRatio={'ratio':numIntens/denIntens,'desc':desc,
                'temperature':outTemperature,'eDensity':outDensity,'filename':intensityRatioFileName, 'numIdx':num_idx, 'denIdx':den_idx}
        #
        # -------------------------------------------------------------------------------------
        #
    def intensityRatioSave(self,outFile=0):
        '''
        Save the intensity ratio to a file.

        The intensity ratio as a function to temperature and eDensity is saved to an asciii file.

        Descriptive information is included at the top of the file.
        '''
        if not outFile:
            outFile = self.IntensityRatio['filename']
#            if chInteractive:
            print(' saving ratio to filename = %s'%(outFile))
        if hasattr(self, 'IntensityRatio'):
            temperature=self.IntensityRatio['temperature']
            eDensity=self.IntensityRatio['eDensity']
            ratio=self.IntensityRatio['ratio']
            out = open(outFile,'w')
            nvalues=len(ratio)
            #
            #  need to add 7 lines to maintain IDL like files
            #
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
#            if chInteractive:
            print(' in .intensityRatioSave(), no IntensityRatio is found')
