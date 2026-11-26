"""
FreeFree module

"""

#import os
#
import numpy as np
#from scipy.interpolate import splev, splrep
from scipy.interpolate import interpn
#from scipy import special
#from scipy.ndimage import map_coordinates
##import matplotlib.pyplot as plt
##import ChiantiPy.core as ch
#from .Ion import ion as chion
from ChiantiPy.base import ionTrails
#from ChiantiPy.base import ionTrails
import ChiantiPy.tools.data as chdata
import ChiantiPy.tools.util as util
import ChiantiPy.tools.io as io
import ChiantiPy.tools.constants as const


class freeFree(ionTrails):
    def __init__(self, Z, abundance=None):
        """

        :param Z: nuclear charge
        :type Z: int
        includeAbund : `bool`, optional
            If True, include the ion abundance in the final output.
        em : array-like, optional
        Line-of-sight emission measure (:math:`\\int\\mathrm{d}l\\,n_en_H`), in units of
            :math:`\\mathrm{cm}^{-5}`, or the volumetric emission measure (:math:`\\int\\mathrm{d}V\\,n_en_H`)
        in units of :math:`\\mathrm{cm}^{-3}`.

        References
        ----------
        .. [110] van Hoof, P. A. M., et al., 2015,  MNRAS, 449, 2112 _

        """
        self.Z = Z
        self.Labels = util.units(chdata.Defaults)

        if not hasattr(self, 'VanHoof'):
            self.VanHoof = io.gff_vanHoofRead(Z)

        if abundance is not None:
            if isinstance(abundance,  float):
                self.Abundance = abundance
                # need the following for p2eRatio but will use default abundance
                self.AbundAll = chdata.AbundanceDefault['abundance']
            elif isinstance(abundance, str):
                ab = io.abundanceRead(abundance,  verbose=True)
                self.Abundance = ab['abundance'][self.Z-1]
                self.AbundanceName = abundance
                self.AbundAll = ab['abundance']

        else:
            self.Abundance = chdata.AbundanceDefault['abundance'][self.Z-1]
            self.AbundAll = chdata.AbundanceDefault['abundance']


        self.Defaults = chdata.Defaults

        if self.Defaults['wavelength'] != 'angstrom':
            print(' the continuum can only be calculated for wavelengths in angstroms')
            return

        return

    def freeFree(self, temperature,  wavelength, includeAbund=True,  em=None):
        """

        :param temperature: temperature in K
        :type temperature: float, array-type

        :param wavelength: wavelength in Angstroms
        :type wavelength: float, array-type

        """

        self.argCheck(temperature=temperature, eDensity=None, pDensity=None, em=em)


        xlabel = self.Labels['xlabel']
        ylabel = self.Labels['spectrumYlabel']

        temperature = np.atleast_1d(temperature)
        self.Temperature = temperature
        wavelength = np.atleast_1d(wavelength)
        self.Wavelength = wavelength

        if not hasattr(self, 'VanHoof'):
            self.VanHoof = io.gff_vanHoofRead(self.Z)
        gff = self.VanHoof['gff']
        logU = self.VanHoof['logU']
        logGam2 = self.VanHoof['logGam2']

#        myGam2 = self.Z**2*const.rydbergErg/(const.boltzmann*temperature)
#        myLogGam2 = np.log10(myGam2)

        self.Gff = np.zeros((temperature.shape[0], wavelength.shape[0]), np.float64)

        if self.Ntemp == 1:
            myGam2 = self.Z**2*const.rydbergErg/(const.boltzmann*temperature)
            myLogGam2 = np.log10(myGam2)
            myU = const.planck*(1.e8*const.light)/(const.boltzmann*temperature*wavelength)
            myLogU = np.log10(myU)
            self.Gff = interpn((logU, logGam2), gff, (myLogU, myLogGam2), method = "cubic", bounds_error = False, fill_value = 1.)
        else:
            self.Gff = np.zeros((temperature.shape[0], wavelength.shape[0]), np.float64)
            for itemp, atemp in enumerate(temperature):
                myGam2 = self.Z**2*const.rydbergErg/(const.boltzmann*atemp)
                myLogGam2 = np.log10(myGam2)
                myU = const.planck*(1.e8*const.light)/(const.boltzmann*atemp*wavelength)
                myLogU = np.log10(myU)
                self.Gff[itemp] = interpn((logU, logGam2), gff, (myLogU, myLogGam2), method = "cubic", bounds_error = False,
                    fill_value = 1.)

        # define the numerical prefactor
        prefactor = ((const.light*1e8)/3./const.emass
                     * (const.fine*const.planck/np.pi)**3
                     * np.sqrt(2.*np.pi/3./const.emass/const.boltzmann))
        # include temperature dependence

        prefactor *= self.Z**2/np.sqrt(temperature)
        if includeAbund:
            prefactor *= self.Abundance
        if self.Em is not None:
            prefactor *= self.Em

        # express in units of ergs or photons
        energy_factor = 1.0
        if chdata.Defaults['flux'] == 'photon':
            energy_factor = const.planck*(1.e8*const.light)/wavelength
        # define exponential factor
        exp_factor = np.exp(-const.planck*(1.e8*const.light)/const.boltzmann
                            / np.outer(self.Temperature, wavelength))/(wavelength**2)

        if prefactor.size > 1:
            free_free_emission = (prefactor[:,np.newaxis]*exp_factor*self.Gff/energy_factor)
        else:
            free_free_emission = (prefactor*exp_factor*self.Gff/energy_factor).squeeze()
        self.FreeFree = {'intensity':free_free_emission, 'temperature':self.Temperature, 'wvl':wavelength,
            'em':self.Em, 'Z':self.Z, 'xlabel':xlabel,  'ylabel':ylabel, 'prefactor':prefactor,
            'exp_factor':exp_factor}



