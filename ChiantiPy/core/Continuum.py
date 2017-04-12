"""
Continuum module
"""
import os

import numpy as np
from scipy.interpolate import splev, splrep
from scipy.ndimage import map_coordinates

import ChiantiPy.tools.data as ch_data
import ChiantiPy.tools.util as ch_util
import ChiantiPy.tools.io as ch_io
import ChiantiPy.tools.constants as ch_const
import ChiantiPy.Gui as chGui


class Continuum:
    """
    The top level class for continuum calculations. Includes methods for the calculation of the
    free-free and free-bound continua.

    Parameters
    ----------
    ionStr : `str`
        CHIANTI notation for the given ion, e.g. 'fe_12' that corresponds to the `Fe XII` ion.
    temperature : array-like
        Temperature  (Kelvin)
    abundance : `float` or `str`, optional
        Elemental abundance relative to Hydrogen or name of CHIANTI abundance file,
        without the '.abund' suffix, e.g. 'sun_photospheric_1998_grevesse'.
    emission_measure : array-like, optional
        Line-of-sight emission measure (:math:`\int\mathrm{d}l\,n_en_H`), in units of
        :math:`\mathrm{cm}^{-5}`, or the volumetric emission measure (:math:`\int\mathrm{d}V\,n_en_H`)
        in units of :math:`\mathrm{cm}^{-3}`.
    """

    def __init__(self, ion_string,  temperature, abundance=None, emission_measure=None):
        nameDict = ch_util.convertName(ion_string)
        self.Z = nameDict['Z']
        self.z_ion = nameDict['Ion']
        self.temperature = np.array(temperature)

        # Throw exception if neutral
        if self.z_ion == 1:
            raise ValueError('{} is a neutral ion and does not produce a continuum'.format(ion_string))

        # Read ionization potential
        self.ionization_potential = ch_data.Ip[self.Z-1, self.z_ion-1]*ch_const.ev2Erg
        # Set abundance
        if abundance is not None:
            try:
                self.abundance = float(abundance)
            except ValueError:
                if abundance in ch_data.AbundanceList:
                    self.abundance_name = abundance
                else:
                    abundChoices = ch_data.AbundanceList
                    abundChoice = chGui.gui.selectorDialog(abundChoices, label='Select Abundance name')
                    abundChoice_idx = abundChoice.selectedIndex
                    self.abundance_name = abundChoices[abundChoice_idx[0]]
        else:
            self.abundance_name = ch_data.Defaults['abundfile']
        if not hasattr(self, 'abundance'):
            self.abundance = ch_data.Abundance[self.abundance_name]['abundance'][self.Z-1]

        # Make sure emission measure is an array if not None
        if emission_measure:
            self.emission_measure = np.array(emission_measure)
        else:
            self.emission_measure = emission_measure

    def free_free(self, wavelength, include_abundance=True, include_ioneq=True, **kwargs):
        """
        Calculates the free-free emission for a single ion.

        Parameters
        ----------
        wavelength : array-like
            In units of angstroms
        include_abundance : `bool`, optional
            If True, include the ion abundance in the final output.
        include_ioneq : `bool`, optional
            If True, include the ionization equilibrium in the final output

        The free-free emission for the given ion is calculated according Eq. [some no] of [1]_,

                    Need some equations here...

        where .... The gaunt factors, , are estimated using the methods of [2]_ and [3]_,
        depending on the temperature and energy regime.

        Notes
        -----
        Can include elemental abundance and ionization equilibrium population
        and the emission measure if specified.
        """
        # define the numerical prefactor
        prefactor = (1.e8*ch_const.light/3./ch_const.emass
                     * np.sqrt(2.*np.pi/3./ch_const.emass/ch_const.boltzmann)
                     * (ch_const.fine*ch_const.planck/np.pi)**3)
        # include temperature dependence
        prefactor *= self.Z**2/np.sqrt(self.temperature)
        if include_abundance:
            prefactor *= self.abundance
        if include_ioneq:
            prefactor *= self.ioneq_one(**kwargs)
        # define exponential factor
        exp_factor = np.exp(-ch_const.planck*(1.e8*ch_const.light)/ch_const.boltzmann
                            / np.outer(self.temperature, wavelength))/(wavelength**2)
        # calculate gaunt factor
        gf_itoh = self.itoh_gaunt_factor(wavelength)
        gf_sutherland = self.sutherland_gaunt_factor(wavelength)
        gf = np.where(np.isnan(gf_itoh), gf_sutherland, gf_itoh)
        # express in units of ergs or photons
        energy_factor = 1.0
        if ch_data.Defaults['flux'] == 'photon':
            energy_factor = ch_const.planck*(1.e8*ch_const.light)/wavelength

        return prefactor[:,np.newaxis]*exp_factor*gf/energy_factor

    def itoh_gaunt_factor(self, wavelength):
        """
        Calculates the free-free gaunt factors of [1]_.

            Need some equations here...

        Notes
        -----
        The relativistic values are valid for :math:`6<\log_{10}(T)< 8.5` and :math:`-4<\log_{10}(u)<1`

        References
        ----------
        .. [1] Itoh, N. et al., 2000, ApJS, `128, 125
            <http://adsabs.harvard.edu/abs/2000ApJS..128..125I>`_
        """
        # calculate scaled energy and temperature
        lower_u = ch_const.planck*(1.e8*ch_const.light)/ch_const.boltzmann/np.outer(self.temperature,wavelength)
        upper_u = 1./2.5*(np.log10(lower_u) + 1.5)
        t = 1./1.25*(np.log10(self.temperature) - 7.25)
        # read in Itoh coefficients
        itoh_coefficients = ch_io.itohRead()['itohCoef'][self.Z - 1].reshape(11,11)
        # calculate Gaunt factor
        gf = np.zeros(upper_u.shape)
        for j in range(11):
            for i in range(11):
                gf += (itoh_coefficients[i,j]*(t**i))[:,np.newaxis]*(upper_u**j)
        # apply NaNs where Itoh approximation is not valid
        gf = np.where(np.logical_and(np.log10(lower_u) >= -4., np.log10(lower_u) <= 1.0),gf,np.nan)
        gf[np.where(np.logical_or(np.log10(self.temperature) <= 6.0,
                                  np.log10(self.temperature) >= 8.5)),:] = np.nan

        return gf

    def sutherland_gaunt_factor(self, wavelength):
        """
        Calculates the free-free gaunt factor calculations of [1]_.

            Need some equations here.

        References
        ----------
        .. [1] Sutherland, R. S., 1998, MNRAS, `300, 321
            <http://adsabs.harvard.edu/abs/1998MNRAS.300..321S>`_
        """
        # calculate scaled quantities
        lower_u = ch_const.planck*(1.e8*ch_const.light)/ch_const.boltzmann/np.outer(self.temperature,wavelength)
        gamma_squared = (self.Z**2)*ch_const.ryd2erg/ch_const.boltzmann/self.temperature[:,np.newaxis]*np.ones(lower_u.shape)
        # convert to index coordinates
        i_lower_u = (np.log10(lower_u) + 4.)*10.
        i_gamma_squared = (np.log10(gamma_squared) + 4.)*5.
        # read in sutherland data
        gf_sutherland_data = ch_io.gffRead()
        # interpolate data to scaled quantities
        gf_sutherland = map_coordinates(gf_sutherland_data['gff'],
                                        [i_gamma_squared.flatten(), i_lower_u.flatten()]).reshape(lower_u.shape)

        return np.where(gf_sutherland < 0., 0., gf_sutherland)

    def free_bound(self, wavelength, include_abundance=True, include_ioneq=True, use_verner=True, **kwargs):
        """
        Calculates the free-bound (radiative recombination) continuum emissivity of an ion.
        Provides emissivity in units of ergs :math:`\mathrm{cm}^{-2}` :math:`\mathrm{s}^{-1}`
        :math:`\mathrm{str}^{-1}` :math:`\mathrm{\AA}^{-1}` for an individual ion.

        Need some equations here.

        Notes
        -----
        - Uses the Gaunt factors of [1]_ for recombination to the ground level
        - Uses the photoionization cross sections of [2]_ to develop the free-bound cross section
        - Does not include the elemental abundance or ionization fraction
        - The specified ion is the target ion

        References
        ----------
        .. [1] Karzas and Latter, 1961, ApJSS, `6, 167
            <http://adsabs.harvard.edu/abs/1961ApJS....6..167K>`_
        .. [2] Verner & Yakovlev, 1995, A&AS, `109, 125
            <http://adsabs.harvard.edu/abs/1995A%26AS..109..125V>`_
        """
        # calculate the photon energy in erg
        photon_energy = ch_const.planck*(1.e8*ch_const.light)/wavelength
        prefactor = (2./np.sqrt(2.*np.pi)/(4.*np.pi)/(ch_const.planck*(ch_const.light**3)
                     * (ch_const.emass*ch_const.boltzmann)**(3./2.)))
        # read the free-bound level information for the recombined and recombining ion
        recombined_fblvl = ch_io.fblvlRead('.'.join([ch_util.zion2filename(self.Z, self.z_ion), 'fblvl']))
        recombining_fblvl = ch_io.fblvlRead('.'.join([ch_util.zion2filename(self.Z, self.z_ion+1), 'fblvl']))
        if 'errorMessage' in recombined_fblvl:
            raise ValueError('No free-bound information available for {}'.format(ch_util.zion2name(self.Z, self.z_ion)))
        # get the multiplicity of the ground state of the recombining ion
        if 'errorMessage' in recombining_fblvl:
            omega_0 = 1.
        else:
            omega_0 = recombining_fblvl['mult'][0]

        energy_over_temp_factor = np.outer(1./(self.temperature**1.5), photon_energy**5.)
        # sum over levels of the recombined ion
        sum_factor = np.zeros((len(self.temperature), len(wavelength)))
        for i,omega_i in enumerate(recombined_fblvl['mult']):
            # ionization potential for level i
            ip = self.ionization_potential - recombined_fblvl['ecm'][i]*ch_const.planck*ch_const.light
            # skip level if photon energy is not sufficiently high
            if ip < 0. or np.all(np.max(photon_energy) < (self.ionization_potential - ip)):
                continue
            # calculate cross-section
            if i == 0 and use_verner:
                cross_section = self.verner_cross_section(photon_energy)
            else:
                cross_section = self.karzas_cross_section(photon_energy, ip,
                                                          recombined_fblvl['pqn'][i],
                                                          recombined_fblvl['l'][i])
            scaled_energy = np.exp(-np.outer(1./(ch_const.boltzmann*self.temperature),
                                             photon_energy - ip))
            # the exponential term can go to infinity for low temperatures
            # but if the cross-section is zero this does not matter
            scaled_energy[:,np.where(cross_section == 0.0)] = 1.0
            sum_factor += omega_i/omega_0*scaled_energy*cross_section

        # combine factors
        fb_emiss = prefactor*energy_over_temp_factor*sum_factor
        # include abundance, ionization equilibrium, photon conversion
        if include_abundance:
            fb_emiss *= self.abundance
        if include_ioneq:
            fb_emiss *= self.ioneq_one(**kwargs)[:,np.newaxis]
        if ch_data.Defaults['flux'] == 'photon':
            fb_emiss /= photon_energy

        return fb_emiss

    def verner_cross_section(self, photon_energy):
        """
        Calculates the photoionization cross section using data from [1]_.

        Need some equations here...

        Notes
        -----
        The cross section refers to the next lower ionization stage

        References
        ----------
        .. [1] Verner & Yakovlev, 1995, A&AS, `109, 125
            <http://adsabs.harvard.edu/abs/1995A%26AS..109..125V>`_
        """
        # read verner data
        verner_info = ch_io.vernerRead()
        eth = verner_info['eth'][self.Z,self.z_ion]*ch_const.ev2Erg
        yw = verner_info['yw'][self.Z,self.z_ion]
        ya = verner_info['ya'][self.Z,self.z_ion]
        p = verner_info['p'][self.Z,self.z_ion]
        # convert from megabarn to cm^2
        sigma0 = verner_info['sig0'][self.Z,self.z_ion]*1e-18
        e0 = verner_info['e0'][self.Z,self.z_ion]*ch_const.ev2Erg
        q = 5.5 + verner_info['l'][self.Z,self.z_ion] - 0.5*p

        # scaled photon energy
        y = photon_energy/e0
        # fitting function
        F = ((y - 1.)**2 + yw**2)*(y**(-q))*(1. + np.sqrt(y/ya))**(-p)
        cross_section = sigma0*F

        return np.where(photon_energy < eth, 0., cross_section)

    def karzas_cross_section(self, photon_energy, ionization_potential, n, l):
        """
        Calculate the K&L photoionization cross-sections.

        Need some equations here...
        """
        # numerical constant, in Mbarn
        kl_constant = 1.077294e-1*8065.54e3
        # read in KL gaunt factor data
        karzas_info = ch_io.klgfbRead()
        if n <= karzas_info['klgfb'].shape[0]:
            scaled_energy = np.log10(photon_energy/ionization_potential)
            f_gf = splrep(karzas_info['pe'], karzas_info['klgfb'][n-1,l,:])
            gaunt_factor = np.exp(splev(scaled_energy, f_gf))
        else:
            gaunt_factor = 1.

        # scaled energy factor, converted to cm^-1
        energy_factor = (((ionization_potential/ch_const.planck/ch_const.light)**2.)
                         * ((photon_energy/ch_const.planck/ch_const.light)**(-3)))
        # cross-section, convert to cm^2
        cross_section = (kl_constant*energy_factor*gaunt_factor/n)*1e-18

        return np.where(photon_energy >= ionization_potential, cross_section, 0.)

    def ioneq_one(self, **kwargs):
        '''
        Provide the ionization equilibrium for the selected ion as a function of temperature.
        returned in self.IoneqOne
        this is a duplicate of the method ion.ioneqOne
        '''

        ioneqAll = ch_io.ioneqRead(ioneqname=kwargs.get('ioneqfile', ch_data.Defaults['ioneqfile']))
        ioneqTemperature = ioneqAll['ioneqTemperature']
        ioneq_one = np.zeros_like(self.temperature)
        thisIoneq = ioneqAll['ioneqAll'][self.Z-1,self.z_ion-1].squeeze()

        gioneq = thisIoneq > 0.
        goodt1 = self.temperature >= ioneqTemperature[gioneq].min()
        goodt2 = self.temperature <= ioneqTemperature[gioneq].max()
        goodt = np.logical_and(goodt1, goodt2)
        y2 = splrep(np.log(ioneqTemperature[gioneq]), np.log(thisIoneq[gioneq]), s=0)

        if goodt.sum() > 0:
            if self.temperature.size > 1:
                gIoneq = splev(np.log(self.temperature[goodt]), y2)   #,der=0)
                ioneq_one[goodt] = np.exp(gIoneq)
            else:
                gIoneq = splev(np.log(self.temperature), y2)
                ioneq_one = np.exp(gIoneq)
        else:
            ioneq_one = 0.

        return ioneq_one
