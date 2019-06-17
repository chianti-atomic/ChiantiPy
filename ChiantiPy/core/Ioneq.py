"""
Ionization equilibrium class
"""
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

import ChiantiPy.tools.util as util
import ChiantiPy.tools.io as io
import ChiantiPy.tools.constants as const
import ChiantiPy.tools.data as chdata

from .Ion import ion


class ioneq(object):
    """
    Calculation ion fractions as a function of temperature assuming
    ionization equilibrium.

    Parameters
    ----------
    el_or_z : `int` or `str`
        Atomic number or symbol

    Note
    ----
    When either loading or calculating a set of ion fractions, the temperature and
    ion fractions are returned to the `Temperature` and `Ioneq` attributes, respectively.
    """

    def __init__(self, el_or_z):
        if type(el_or_z) is str:
            self.Z = util.el2z(el_or_z)
        else:
            self.Z = el_or_z

    def load(self, ioneqName = None):
        """
        Read temperature and ion fractions from a CHIANTI ".ioneq" file.
        """
        if ioneqName is None:
            ioneqName = chdata.Defaults['ioneqfile']
        ioneqAll = io.ioneqRead(ioneqName)
        self.Temperature = ioneqAll['ioneqTemperature']
        self.Ioneq = ioneqAll['ioneqAll'][self.Z - 1]
        self.IoneqName = ioneqAll['ioneqname']

    def calculate(self, temperature):
        """
        Calculate ion fractions for given temperature array using the total
        ionization and recombination rates.
        """
        self.Temperature = np.array(temperature, np.float64)
        if self.Temperature.size == 1:
            print(' temperature must be an array')
            return
        ionList = []
        chIons = []
        z = self.Z
        for stage in range(1, z+2):
            ionStr = util.zion2name(z, stage)
            ionList.append(ionStr)
            atom = ion(ionStr, temperature = self.Temperature, setup=0)
            atom.setupIonrec()
            atom.ionizRate()
            atom.recombRate()
            chIons.append(atom)

        ntemp = chIons[0].IonizRate['temperature'].size
        if ntemp == 1:
            ioneq = np.zeros((z+1), np.float64)
            factor = []
            for anIon in chIons:
                if hasattr(anIon, 'RecombRate') and hasattr(anIon, 'IonizRate'):
                    rat = anIon.IonizRate['rate']/anIon.RecombRate['rate']
                    factor.append(rat**2 + rat**(-2))
                else:
                    factor.append(0.)
            factor[0] = max(factor)
            factor[-1] = max(factor)
            ionmax = factor.index(min(factor))
            ioneq[ionmax] = 1.

            for iz in range(ionmax+1, z+1):
                ionrate = chIons[iz-1].IonizRate['rate']
                recrate = chIons[iz].RecombRate['rate']
                ioneq[iz] = ionrate*ioneq[iz-1]/recrate

            for iz in range(ionmax-1, -1, -1):
                ionrate = chIons[iz].IonizRate['rate']
                recrate = chIons[iz+1].RecombRate['rate']
                ioneq[iz] = recrate*ioneq[iz+1]/ionrate
            ionsum = ioneq.sum()
            ioneq = ioneq/ionsum
            self.Ioneq = ioneq
        else:
            ioneq = np.zeros((z+1,ntemp ), np.float64)
            for it in range(ntemp):
                factor = []
                for anIon in chIons:
                    if type(anIon.IonizRate) != type(None) and type(anIon.RecombRate) != type(None):
                        ioniz = anIon.IonizRate['rate'][it]
                        recomb = anIon.RecombRate['rate'][it]
                        if ioniz == 0. or recomb == 0.:
                            rat = 1.e-100
                        else:
                            rat = anIon.IonizRate['rate'][it]/anIon.RecombRate['rate'][it]
                        try:
                            factor.append(rat**2 + rat**(-2))
                        except:
                            factor.append(0.)
                    else:
                        factor.append(0.)
                factor[0] = max(factor)
                factor[-1] = max(factor)
                ionmax = factor.index(min(factor))
                ioneq[ionmax, it] = 1.

                for iz in range(ionmax+1, z+1):
                    ionrate = chIons[iz-1].IonizRate['rate'][it]
                    recrate = chIons[iz].RecombRate['rate'][it]
                    if recrate != 0.:
                        ioneq[iz, it] = ionrate*ioneq[iz-1, it]/recrate
                    else:
                        ioneq[iz, it] = 0.

                for iz in range(ionmax-1, -1, -1):
                    ionrate = chIons[iz].IonizRate['rate'][it]
                    recrate = chIons[iz+1].RecombRate['rate'][it]
                    if ionrate != 0.:
                        ioneq[iz, it] = recrate*ioneq[iz+1, it]/ionrate
                    else:
                        ioneq[iz, it] = 0.
                ionsum = ioneq[:, it].sum()
                ioneq[:, it] = ioneq[:, it]/ionsum
            self.Ioneq = ioneq

    def plot(self, stages=0, xr=0, yr=0, oplot=False, label=1, title=1,  bw=0, semilogx = 0, verbose=0):
        '''
        Plots the ionization equilibria.

        self.plot(xr=None, yr=None, oplot=False)
        stages = sequence of ions to be plotted, neutral == 1, fully stripped == Z+1
        xr = temperature range, yr = ion fraction range

        for overplotting:
        oplot="ioneqfilename" such as 'mazzotta'
        or if oplot=True or oplot=1 and a widget will come up so that a file can be selected.
        '''
        if hasattr(self, 'Ioneq'):
            ioneq = getattr(self, 'Ioneq')
        else:
            print(' must first load or calculate and ionization equilibrium')
            return

        if bw:
            linestyle = ['k-','k--', 'k-.', 'k:']
            plt.rcParams['font.size'] = 16.
            lw = 2
        else:
            linestyle = ['b-','r--', 'g-.', 'm:']
            plt.rcParams['font.size'] = 14.
            lw = 2
        #
        if not stages:
            stages = range(1, self.Z+2)
        elif min(stages) < 1 or max(stages) > self.Z+1:
            stages = range(1, self.Z+2)  #  spectroscopic notation
        if not xr:
            xr = [self.Temperature.min(), self.Temperature.max()]
        if not yr:
            yr = [0.01, 1.1]
        xyr = list(xr)
        xyr.extend(list(yr))
        #
        iz = stages[0]
        if semilogx:
            plt.semilogx(self.Temperature, ioneq[iz-1], linestyle[0], lw=lw)
        else:
            plt.loglog(self.Temperature, ioneq[iz-1], linestyle[0], lw=lw)
        if label:
            idx = self.Ioneq[iz-1] == self.Ioneq[iz-1].max()
            if idx.sum() > 1:
                jdx = np.arange(len(idx))
                idx = int(jdx[idx].max())
            ann = const.Ionstage[iz-1]
            plt.annotate(ann, [self.Temperature[idx], 0.7*ioneq[iz-1, idx]], ha='center')
        for iz in stages[1:]:
            if semilogx:
                plt.semilogx(self.Temperature, ioneq[iz-1], linestyle[0], lw=lw)
            else:
                plt.loglog(self.Temperature, ioneq[iz-1], linestyle[0], lw=lw)
            if label:
                idx = ioneq[iz-1] == ioneq[iz-1].max()
                if idx.sum() > 1:
                    jdx = np.arange(len(idx))
                    idx = int(jdx[idx].mean())
                ann = const.Ionstage[iz-1]
                plt.annotate(ann, [self.Temperature[idx], 0.7*ioneq[iz-1, idx]], ha='center')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Ion Fraction')
        atitle = 'Chianti Ionization Equilibrium for '+const.El[self.Z-1].capitalize()
        #
        if oplot:
            if isinstance(oplot,int):
                result = io.ioneqRead(ioneqName='')
#                print('keys = ', result.keys()
                if result != False:
                    atitle += '  & '+ result['ioneqname'].replace('.ioneq', '')
                    atitle += ' '+linestyle[0]
                    for iz in stages:
                        plt.plot(result['ioneqTemperature'], result['ioneqAll'][self.Z-1, iz-1],linestyle[1], lw=lw)
            elif type(oplot) == type('string'):
                atitle += '  & ' + oplot
                result = io.ioneqRead(ioneqname=oplot)
                if result != False:
                    for iz in stages:
                        plt.plot(result['ioneqTemperature'], result['ioneqAll'][self.Z-1, iz-1],linestyle[1], lw=lw)
            elif type(oplot) == type([]):
                for iplot in range(len(oplot)):
                    result = io.ioneqRead(ioneqname=oplot[iplot])
#                    print 'keys = ', result.keys()
                    if result != False:
                        atitle += '  & '+oplot[iplot]+' '+linestyle[iplot%3]
                        for iz in stages:
                            plt.plot(result['ioneqTemperature'], result['ioneqAll'][self.Z-1, iz-1],linestyle[1], lw=lw)
            else:
                print(' oplot file not understood %s'%(oplot))
        if title:
            plt.title(atitle)
        plt.tight_layout()
        plt.axis(xyr)

    def plotRatio(self, stageN, stageD, xr=0, yr=0, label=1, title=1,  bw=0, semilogx = 1, verbose=0):
        '''
        Plots the ratio of the ionization equilibria of two stages of a given element

        self.plotRatio(stageN, stageD)
        stages = sequence of ions to be plotted, neutral == 1, fully stripped == Z+1
        stageN = numerator
        stageD = denominator
        xr = temperature range, yr = ion fraction range

        '''
        ionN = util.zion2name(self.Z, stageN)
        ionD = util.zion2name(self.Z, stageD)
        ionNS = util.zion2spectroscopic(self.Z, stageN)
        ionDS = util.zion2spectroscopic(self.Z, stageD)

        atitle = 'ratio of %s to %s '%(ionN, ionD)
        alabel = 'ratio of %s to %s '%(ionNS, ionDS)
        print(atitle)
        if not hasattr(self, 'Ioneq'):
            print(' must first load or calculate an ionization equilibrium')
            return

        if bw:
            linestyle = ['k-','k--', 'k-.', 'k:']
            plt.rcParams['font.size'] = 16.
            lw = 2
        else:
            linestyle = ['b-','r--', 'g-.', 'm:']
            plt.rcParams['font.size'] = 14.
            lw = 2
        #
        goodTn = self.Ioneq[stageN - 1, :] > 0.
        goodTd = self.Ioneq[stageD -1, : ] > 0.
        realGoodT = np.logical_and(goodTn, goodTd)
        goodT = self.Temperature[realGoodT]
        goodR = self.Ioneq[stageN - 1, realGoodT]/self.Ioneq[stageD - 1, realGoodT]
        if not xr:
            xr = [goodT.min(), goodT.max()]
        if not yr:
            yr = [0.01, 1.1]
        xyr = list(xr)
        xyr.extend(list(yr))
        #
        if semilogx:
            plt.semilogx(goodT, goodR, linestyle[0], lw=lw, label=alabel)
        else:
            plt.loglog(goodT, goodR, linestyle[0], lw=lw, label=alabel)
        plt.xlabel('Temperature (K)')
        plt.ylabel('Ratio')
        if title:
            atitle = 'CHIANTI Ionization Equilibrium'
            plt.title(atitle)
        plt.legend(loc='lower right')
        plt.tight_layout()
        self.Ratio={'Temperature':goodT, 'Ratio':goodR, 'label':alabel}
#        plt.axis(xyr)
