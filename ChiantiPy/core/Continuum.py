"""
Continuum module
"""
import os

import numpy as np
from scipy.interpolate import splev, splrep
from scipy.ndimage import map_coordinates

from .Ioneq import ioneq
from ChiantiPy.base import ionTrails
import ChiantiPy.tools.data as chdata
import ChiantiPy.tools.util as util
import ChiantiPy.tools.io as io
import ChiantiPy.tools.constants as const
import ChiantiPy.Gui as chGui


class continuum(ionTrails):
    """
    The top level class for continuum calculations. Includes methods for the calculation of the
    free-free and free-bound continua.

    Parameters
    ----------
    ionStr : `str`
        CHIANTI notation for the given ion, e.g. 'fe_12' that corresponds to the Fe XII ion.
    temperature : array-like
        In units of Kelvin
    abundance : `float` or `str`, optional
        Elemental abundance relative to Hydrogen or name of CHIANTI abundance file,
        without the '.abund' suffix, e.g. 'sun_photospheric_1998_grevesse'.
    em : array-like, optional
        Line-of-sight emission measure (:math:`\int\mathrm{d}l\,n_en_H`), in units of
        :math:`\mathrm{cm}^{-5}`, or the volumetric emission measure (:math:`\int\mathrm{d}V\,n_en_H`)
        in units of :math:`\mathrm{cm}^{-3}`.

    Examples
    --------
    >>> import ChiantiPy.core as ch
    >>> import numpy as np
    >>> temperature = np.logspace(4,9,20)
    >>> cont = ch.continuum('fe_15',temperature)
    >>> wavelength = np.arange(1,1000,10)
    >>> cont.freeFree(wavelength)
    >>> cont.freeBound(wavelength, include_abundance=True, include_ioneq=False)
    >>> cont.calculate_free_free_loss()
    >>> cont.calculate_free_bound_loss()

    Notes
    -----
    The methods for calculating the free-free and free-bound emission and losses return
    their result to an attribute. See the respective docstrings for more information.
    """

    def __init__(self, ionStr,  temperature, abundance=None, em=None, verbose=0):
        self.IonStr = ionStr
        self.nameDict = util.convertName(ionStr)
        self.Z = self.nameDict['Z']
        self.Stage = self.nameDict['Ion']
        self.Ion = self.nameDict['Ion']
        self.Temperature = np.atleast_1d(temperature)
        self.NTemperature = self.Temperature.size
        self.Defaults = chdata.Defaults


        self.argCheck(temperature=temperature, eDensity=None, pDensity=None, em=em, verbose=verbose)

        self.Ip = chdata.Ip[self.Z-1, self.Stage-1]
        self.Ipr = chdata.Ip[self.Z-1, self.Stage-2]
        self.ionization_potential = chdata.Ip[self.Z-1, self.Stage-1]*const.ev2Erg
        self.IprErg = self.Ipr*const.ev2Erg
        # Set abundance
        if abundance is not None:
            try:
                self.Abundance = float(abundance)
            except ValueError:
                if abundance in chdata.AbundanceList:
                    self.AbundanceName = abundance
                else:
                    abundChoices = chdata.AbundanceList
                    abundChoice = chGui.gui.selectorDialog(abundChoices, label='Select Abundance name')
                    abundChoice_idx = abundChoice.selectedIndex
                    self.AbundanceName = abundChoices[abundChoice_idx[0]]
        else:
            self.AbundanceName = chdata.Defaults['abundfile']
        if not hasattr(self, 'Abundance'):
            self.Abundance = chdata.Abundance[self.AbundanceName]['abundance'][self.Z-1]
        self.ioneqOne()

    def free_free_loss(self, **kwargs):
        """
        Calculate the free-free energy loss rate of an ion. The result is returned to the
        `free_free_loss` attribute.

        The free-free radiative loss rate is given by Eq. 5.15a of [1]_. Writing the numerical
        constant in terms of the fine structure constant :math:`\\alpha`,

        .. math::
           \\frac{dW}{dtdV} = \\frac{4\\alpha^3h^2}{3\pi^2m_e}\left(\\frac{2\pi k_B}{3m_e}\\right)^{1/2}Z^2T^{1/2}\\bar{g}_B

        where where :math:`Z` is the nuclear charge, :math:`T` is the electron temperature, and
        :math:`\\bar{g}_{B}` is the wavelength-averaged and velocity-averaged Gaunt factor. The
        Gaunt factor is calculated using the methods of [2]_. Note that this expression for the
        loss rate is just the integral over wavelength of Eq. 5.14a of [1]_, the free-free emission, and
        is expressed in units of erg :math:`\mathrm{cm}^3\,\mathrm{s}^{-1}`.

        References
        ----------
        .. [1] Rybicki and Lightman, 1979, Radiative Processes in Astrophysics,
            `(Wiley-VCH) <http://adsabs.harvard.edu/abs/1986rpa..book.....R>`_
        .. [2] Karzas and Latter, 1961, ApJSS, `6, 167
            <http://adsabs.harvard.edu/abs/1961ApJS....6..167K>`_
        """
        # interpolate wavelength-averaged K&L gaunt factors
        gf_kl_info = io.gffintRead()
        gamma_squared = self.ionization_potential/const.boltzmann/self.Temperature
        for i, atemp in enumerate(self.Temperature):
            print('%s T:  %10.2e gamma_squared  %10.2e'%(self.IonStr, atemp, gamma_squared[i]))
        gaunt_factor = splev(np.log(gamma_squared),
                             splrep(gf_kl_info['g2'],gf_kl_info['gffint']), ext=3)
        # calculate numerical constant
        prefactor = (4.*(const.fine**3)*(const.planck**2)/3./(np.pi**2)/const.emass
                     * np.sqrt(2.*np.pi*const.boltzmann/3./const.emass))
        # include abundance and ionization equilibrium
        prefactor *= self.Abundance*self.ioneq_one(self.Stage, **kwargs)

        self.free_free_loss = prefactor*(self.Z**2)*np.sqrt(self.Temperature)*gaunt_factor

    def freeFree(self, wavelength, include_abundance=True, include_ioneq=True, **kwargs):
        """
        Calculates the free-free emission for a single ion. The result is returned as a dict to
        the `FreeFree` attribute.  The dict has the keywords `intensity`, `wvl`, `temperature`, `em`.

        The free-free emission for the given ion is calculated according Eq. 5.14a of [1]_,
        substituting :math:`\\nu=c/\lambda`, dividing by the solid angle, and writing the numerical
        constant in terms of the fine structure constant :math:`\\alpha`,

        .. math::
           \\frac{dW}{dtdVd\lambda} = \\frac{c}{3m_e}\left(\\frac{\\alpha h}{\pi}\\right)^3\left(\\frac{2\pi}{3m_ek_B}\\right)^{1/2}\\frac{Z^2}{\lambda^2T^{1/2}}\exp{\left(-\\frac{hc}{\lambda k_BT}\\right)}\\bar{g}_{ff},

        where :math:`Z` is the nuclear charge, :math:`T` is the electron temperature in K, and
        :math:`\\bar{g}_{ff}` is the velocity-averaged Gaunt factor. The Gaunt factor is estimated
        using the methods of [2]_ and [3]_, depending on the temperature and energy regime. See
        `itoh_gaunt_factor` and `sutherland_gaunt_factor` for more details.

        The free-free emission is in units of erg
        :math:`\mathrm{cm}^3\mathrm{s}^{-1}\mathrm{\mathring{A}}^{-1}\mathrm{str}^{-1}`. If the emission
        measure has been set, the units will be multiplied by :math:`\mathrm{cm}^{-5}` or
        :math:`\mathrm{cm}^{-3}`, depending on whether it is the line-of-sight or volumetric
        emission measure, respectively.

        Parameters
        ----------
        wavelength : array-like
            In units of angstroms
        include_abundance : `bool`, optional
            If True, include the ion abundance in the final output.
        include_ioneq : `bool`, optional
            If True, include the ionization equilibrium in the final output

        References
        ----------
        .. [1] Rybicki and Lightman, 1979, Radiative Processes in Astrophysics,
            `(Wiley-VCH) <http://adsabs.harvard.edu/abs/1986rpa..book.....R>`_
        .. [2] Itoh, N. et al., 2000, ApJS, `128, 125
            <http://adsabs.harvard.edu/abs/2000ApJS..128..125I>`_
        .. [3] Sutherland, R. S., 1998, MNRAS, `300, 321
            <http://adsabs.harvard.edu/abs/1998MNRAS.300..321S>`_
        """
        wavelength = np.atleast_1d(wavelength)
        # define the numerical prefactor
        prefactor = ((const.light*1e8)/3./const.emass
                     * (const.fine*const.planck/np.pi)**3
                     * np.sqrt(2.*np.pi/3./const.emass/const.boltzmann))
        # include temperature dependence
        prefactor *= self.Z**2/np.sqrt(self.Temperature)
        if include_abundance:
            prefactor *= self.Abundance
        if include_ioneq:
            prefactor *= self.IoneqOne
        if self.Em is not None:
            prefactor *= self.Em
        # define exponential factor
        exp_factor = np.exp(-const.planck*(1.e8*const.light)/const.boltzmann
                            / np.outer(self.Temperature, wavelength))/(wavelength**2)
        # calculate gaunt factor
        gf_itoh = self.itoh_gaunt_factor(wavelength)
        gf_sutherland = self.sutherland_gaunt_factor(wavelength)
        gf = np.where(np.isnan(gf_itoh), gf_sutherland, gf_itoh)
        # express in units of ergs or photons
        energy_factor = 1.0
        if chdata.Defaults['flux'] == 'photon':
            energy_factor = const.planck*(1.e8*const.light)/wavelength

        free_free_emission = (prefactor[:,np.newaxis]*exp_factor*gf/energy_factor).squeeze()
        self.FreeFree = {'intensity':free_free_emission, 'temperature':self.Temperature, 'wvl':wavelength, 'em':self.Em, 'ions':self.IonStr}

    def freeFreeLoss(self, **kwargs):
        """
        Calculate the free-free energy loss rate of an ion. The result is returned to the
        `FreeFreeLoss` attribute.

        The free-free radiative loss rate is given by Eq. 5.15a of [1]_. Writing the numerical
        constant in terms of the fine structure constant :math:`\\alpha`,

        .. math::
           \\frac{dW}{dtdV} = \\frac{4\\alpha^3h^2}{3\pi^2m_e}\left(\\frac{2\pi k_B}{3m_e}\\right)^{1/2}Z^2T^{1/2}\\bar{g}_B

        where where :math:`Z` is the nuclear charge, :math:`T` is the electron temperature, and
        :math:`\\bar{g}_{B}` is the wavelength-averaged and velocity-averaged Gaunt factor. The
        Gaunt factor is calculated using the methods of [2]_. Note that this expression for the
        loss rate is just the integral over wavelength of Eq. 5.14a of [1]_, the free-free emission, and
        is expressed in units of erg :math:`\mathrm{cm}^3\,\mathrm{s}^{-1}`.

        References
        ----------
        .. [1] Rybicki and Lightman, 1979, Radiative Processes in Astrophysics,
            `(Wiley-VCH) <http://adsabs.harvard.edu/abs/1986rpa..book.....R>`_
        .. [2] Karzas and Latter, 1961, ApJSS, `6, 167
            <http://adsabs.harvard.edu/abs/1961ApJS....6..167K>`_
        """
        # interpolate wavelength-averaged K&L gaunt factors
        gf_kl_info = io.gffintRead()
        gamma_squared = self.IprErg/const.boltzmann/self.Temperature
#        for i, atemp in enumerate(self.Temperature):
#            print('%s T:  %10.2e gamma_squared  %10.2e'%(self.IonStr, atemp, gamma_squared[i]))
        gaunt_factor = splev(np.log(gamma_squared),
                             splrep(gf_kl_info['g2'],gf_kl_info['gffint']), ext=3)
        # calculate numerical constant
        prefactor = (4.*(const.fine**3)*(const.planck**2)/3./(np.pi**2)/const.emass
                     * np.sqrt(2.*np.pi*const.boltzmann/3./const.emass))
        # include abundance and ionization equilibrium
        prefactor *= self.Abundance*self.IoneqOne

        self.FreeFreeLoss = {'rate':prefactor*(self.Z**2)*np.sqrt(self.Temperature)*gaunt_factor}


    def itoh_gaunt_factor(self, wavelength):
        """
        Calculates the free-free gaunt factors of [1]_.

        An analytic fitting formulae for the relativistic Gaunt factor is given by Eq. 4 of [1]_,

        .. math::
           g_{Z} = \sum^{10}_{i,j=0}a_{ij}t^iU^j

        where,

        .. math::
           t = \\frac{1}{1.25}(\log_{10}{T} - 7.25),\\
           U = \\frac{1}{2.5}(\log_{10}{u} + 1.5),

        :math:`u=hc/\lambda k_BT`, and :math:`a_{ij}` are the fitting coefficients and are read
        in using `ChiantiPy.tools.io.itohRead` and are given in Table 4 of [1]_. These values
        are valid for :math:`6<\log_{10}(T)< 8.5` and :math:`-4<\log_{10}(u)<1`.

        See Also
        --------
        ChiantiPy.tools.io.itohRead : Read in Gaunt factor coefficients from [1]_

        References
        ----------
        .. [1] Itoh, N. et al., 2000, ApJS, `128, 125
            <http://adsabs.harvard.edu/abs/2000ApJS..128..125I>`_
        """
        # calculate scaled energy and temperature
        lower_u = const.planck*(1.e8*const.light)/const.boltzmann/np.outer(self.Temperature, wavelength)
        upper_u = 1./2.5*(np.log10(lower_u) + 1.5)
        t = 1./1.25*(np.log10(self.Temperature) - 7.25)
        # read in Itoh coefficients
        itoh_coefficients = io.itohRead()['itohCoef'][self.Z - 1].reshape(11,11)
        # calculate Gaunt factor
        gf = np.zeros(upper_u.shape)
        for j in range(11):
            for i in range(11):
                gf += (itoh_coefficients[i,j]*(t**i))[:,np.newaxis]*(upper_u**j)
        # apply NaNs where Itoh approximation is not valid
        gf = np.where(np.logical_and(np.log10(lower_u) >= -4., np.log10(lower_u) <= 1.0),gf,np.nan)
        gf[np.where(np.logical_or(np.log10(self.Temperature) <= 6.0,
                                  np.log10(self.Temperature) >= 8.5)),:] = np.nan

        return gf

    def sutherland_gaunt_factor(self, wavelength):
        """
        Calculates the free-free gaunt factor calculations of [1]_.

        The Gaunt factors of [1]_ are read in using `ChiantiPy.tools.io.gffRead`
        as a function of :math:`u` and :math:`\gamma^2`. The data are interpolated
        to the appropriate wavelength and temperature values using
        `~scipy.ndimage.map_coordinates`.

        References
        ----------
        .. [1] Sutherland, R. S., 1998, MNRAS, `300, 321 <http://adsabs.harvard.edu/abs/1998MNRAS.300..321S>`_
        """
        # calculate scaled quantities
        lower_u = const.planck*(1.e8*const.light)/const.boltzmann/np.outer(self.Temperature,wavelength)
        gamma_squared = (self.Z**2)*const.ryd2erg/const.boltzmann/self.Temperature[:,np.newaxis]*np.ones(lower_u.shape)
        # convert to index coordinates
        i_lower_u = (np.log10(lower_u) + 4.)*10.
        i_gamma_squared = (np.log10(gamma_squared) + 4.)*5.
        # read in sutherland data
        gf_sutherland_data = io.gffRead()
        # interpolate data to scaled quantities
        gf_sutherland = map_coordinates(gf_sutherland_data['gff'],
                                        [i_gamma_squared.flatten(), i_lower_u.flatten()]).reshape(lower_u.shape)

        return np.where(gf_sutherland < 0., 0., gf_sutherland)

    def calculate_free_bound_loss(self, **kwargs):
        """
        Calculate the free-bound energy loss rate of an ion. The result is returned to the
        `free_bound_loss` attribute.

        The free-bound loss rate can be calculated by integrating the free-bound emission over the wavelength.
        This is difficult using the expression in `calculate_free_bound_emission` so we instead use the
        approach of [1]_ and [2]_. Eq. 1a of [2]_ can be integrated over wavelength to get the free-bound loss rate,

        .. math::
           \\frac{dW}{dtdV} = C_{ff}\\frac{k}{hc}T^{1/2}G_{fb},

        in units of erg :math:`\mathrm{cm}^3\,\mathrm{s}^{-1}` where :math:`G_{fb}` is the free-bound Gaunt factor as
        given by Eq. 15 of [2]_ (see `mewe_gaunt_factor` for more details) and :math:`C_{ff}` is the numerical constant
        as given in Eq. 4 of [1]_ and can be written in terms of the fine structure constant :math:`\\alpha`,

        .. math::
           C_{ff}\\frac{k}{hc} = \\frac{8}{3}\left(\\frac{\pi}{6}\\right)^{1/2}\\frac{h^2\\alpha^3}{\pi^2}\\frac{k_B}{m_e^{3/2}} \\approx 1.43\\times10^{-27}

        References
        ----------
        .. [1] Gronenschild, E.H.B.M. and Mewe, R., 1978, A&AS, `32, 283 <http://adsabs.harvard.edu/abs/1978A%26AS...32..283G>`_
        .. [2] Mewe, R. et al., 1986, A&AS, `65, 511 <http://adsabs.harvard.edu/abs/1986A%26AS...65..511M>`_
        """
        # Calculate Gaunt factor according to Mewe
        gaunt_factor = self.mewe_gaunt_factor()
        # Numerical prefactor
        prefactor = (8./3.*np.sqrt(np.pi/6.)*(const.planck**2)*(const.fine**3)/(np.pi**2)
                     * (const.boltzmann**(1./2.))/(const.emass**(3./2.)))

        self.free_bound_loss = gaunt_factor*np.sqrt(self.Temperature)*prefactor

    def freeBoundLossMewe(self, **kwargs):
        """
        Calculate the free-bound energy loss rate of an ion. The result is returned to the
        `free_bound_loss` attribute.

        The free-bound loss rate can be calculated by integrating the free-bound emission over the wavelength.
        This is difficult using the expression in `calculate_free_bound_emission` so we instead use the
        approach of [1]_ and [2]_. Eq. 1a of [2]_ can be integrated over wavelength to get the free-bound loss rate,

        .. math::
           \\frac{dW}{dtdV} = C_{ff}\\frac{k}{hc}T^{1/2}G_{fb},

        in units of erg :math:`\mathrm{cm}^3\,\mathrm{s}^{-1}` where :math:`G_{fb}` is the free-bound Gaunt factor as
        given by Eq. 15 of [2]_ (see `mewe_gaunt_factor` for more details) and :math:`C_{ff}` is the numerical constant
        as given in Eq. 4 of [1]_ and can be written in terms of the fine structure constant :math:`\\alpha`,

        .. math::
           C_{ff}\\frac{k}{hc} = \\frac{8}{3}\left(\\frac{\pi}{6}\\right)^{1/2}\\frac{h^2\\alpha^3}{\pi^2}\\frac{k_B}{m_e^{3/2}} \\approx 1.43\\times10^{-27}

        References
        ----------
        .. [1] Gronenschild, E.H.B.M. and Mewe, R., 1978, A&AS, `32, 283 <http://adsabs.harvard.edu/abs/1978A%26AS...32..283G>`_
        .. [2] Mewe, R. et al., 1986, A&AS, `65, 511 <http://adsabs.harvard.edu/abs/1986A%26AS...65..511M>`_
        """
        nameDict = util.convertName(self.IonStr)
        lower = nameDict['lower']
        self.Recombined_fblvl = io.fblvlRead(lower)
        if 'errorMessage' in self.Recombined_fblvl:
            errorMessage = 'No free-bound information available for {}'.format(self.IonStr)
            rate = np.zeros_like(self.Temperature)
            self.FreeBoundLoss = {'rate':rate, 'errorMessage':errorMessage}
            return
# Calculate Gaunt factor according to Mewe
        gaunt_factor = self.mewe_gaunt_factor()
        # Numerical prefactor
        prefactor = (8./3.*np.sqrt(np.pi/6.)*(const.planck**2)*(const.fine**3)/(np.pi**2)
                     * (const.boltzmann**(1./2.))/(const.emass**(3./2.)))

        self.FreeBoundLoss = {'rate':gaunt_factor*np.sqrt(self.Temperature)*prefactor, 'temperature':self.Temperature}

    def mewe_gaunt_factor(self, **kwargs):
        """
        Calculate the Gaunt factor according to [1]_ for a single ion :math:`Z_z`.

        Using Eq. 9 of [1]_, the free-bound Gaunt factor for a single ion can be written as,

        .. math::
           G_{fb}^{Z,z} = \\frac{E_H}{k_BT}\\mathrm{Ab}(Z)\\frac{N(Z,z)}{N(Z)}f(Z,z,n)

        where :math:`E_H` is the ground-state potential of H, :math:`\mathrm{Ab}(Z)` is the
        elemental abundance, :math:`\\frac{N(Z,z)}{N(Z)}` is the fractional ionization, and
        :math:`f(Z,z,n)` is given by Eq. 10 and is approximated by Eq 16 as,

        .. math::
           f(Z,z,n) \\approx f_2(Z,z,n_0) = 0.9\\frac{\zeta_0z_0^4}{n_0^5}\exp{\left(\\frac{E_Hz_0^2}{n_0^2k_BT}\\right)} + 0.42\\frac{z^4}{n_0^{3/2}}\exp{\left(\\frac{E_Hz^2}{(n_0 + 1)^2k_BT}\\right)}

        where :math:`n_0` is the principal quantum number, :math:`z_0` is the effective charge (see Eq. 7 of [1]_),
        and :math:`\zeta_0` is the number of vacancies in the 0th shell and is given in Table 1 of [1]_.
        Here it is calculated in the same manner as in `fb_rad_loss.pro <http://www.chiantidatabase.org/idl/continuum/fb_rad_loss.pro>`_
        of the CHIANTI IDL library. Note that in the expression for :math:`G_{fb}`, we have not included
        the :math:`N_H/n_e` factor.

        Raises
        ------
        ValueError
            If no .fblvl file is available for this ion

        References
        ----------
        .. [1] Mewe, R. et al., 1986, A&AS, `65, 511 <http://adsabs.harvard.edu/abs/1986A%26AS...65..511M>`_
        """
        # read in free-bound level information for the recombined ion
        # thermal energy scaled by H ionization potential
        scaled_energy = const.ryd2erg/const.boltzmann/self.Temperature
        # set variables used in Eq. 16 of Mewe et al.(1986)
        n_0 = self.Recombined_fblvl['pqn'][0]
#        z_0 = np.sqrt(self.ionization_potential/const.ryd2erg)*n_0
        z_0 = np.sqrt(self.Ipr/const.ryd2erg)*n_0

        # calculate zeta_0, the number of vacancies in the recombining ion
        # see zeta_0 function in chianti/idl/continuum/fb_rad_loss.pro and
        # Table 1 of Mewe et al. (1986)
        if self.Z - self.Stage > 22:
            zeta_0 = self.Z - self.Stage + 55
        elif 8 < self.Z - self.Stage <= 22:
            zeta_0 = self.Z - self.Stage + 27
        elif 0 < self.Z - self.Stage <= 8:
            zeta_0 = self.Z - self.Stage + 9
        else:
            zeta_0 = self.Z - self.Stage + 1

        ip = self.Ipr - self.Recombined_fblvl['ecm'][0]*const.planck*const.light
#        ip = self.ionization_potential - recombined_fblvl['ecm'][0]*const.planck*const.light
        f_2 = (0.9*zeta_0*(z_0**4)/(n_0**5)*np.exp(scaled_energy*(z_0**2)/(n_0**2) - ip/const.boltzmann/self.Temperature)
               + 0.42/(n_0**1.5)*(self.Stage**4))

#        return scaled_energy*f_2*self.Abundance*self.ioneq_one(self.Stage+1, **kwargs)
        return scaled_energy*f_2*self.Abundance*self.IoneqOne
            #

    def freeBoundLoss(self):
        '''
        to calculate the free-bound (radiative recombination) energy loss rate coefficient of an ion,
        the ion is taken to be the target ion,
        including the elemental abundance and the ionization equilibrium population
        uses the Gaunt factors of Karzas, W.J, Latter, R, 1961, ApJS, 6, 167
        provides rate = ergs cm^-2 s^-1
        '''
        #
        temperature = self.Temperature
        #
        nameDict = util.convertName(self.IonStr)
        lowerDict = util.convertName(nameDict['lower'])
        if hasattr(self, 'Fblvl'):
            fblvl = self.Fblvl
        else:
            fblvlname = nameDict['filename']+'.fblvl'
            if os.path.isfile(fblvlname):
                self.Fblvl = io.fblvlRead(self.IonStr)
                fblvl = self.Fblvl
            elif self.Stage == self.Z+1:
                fblvl = {'mult':[1., 1.]}
            else:
                self.FreeBoundLoss = {'errorMessage':' file does not exist %s .fblvl'%(fblvlname)}
                return
        #  need some data for the recombined ion
        #
        if hasattr(self, 'rFblvl'):
            rFblvl = self.rFblvl
        else:
            rfblvlname = lowerDict['filename']+'.fblvl'
            if os.path.isfile(rfblvlname):
                self.rFblvl = io.fblvlRead(nameDict['lower'])
                rFblvl = self.rFblvl
            else:
                self.FreeBoundLoss = {'errorMessage':' file does not exist %s .fblvl'%(rfblvlname)}
                return
        #
        gIoneq = self.IoneqOne
        #
        abund = self.Abundance
        #
        #
        nlvls = len(rFblvl['lvl'])
        # pqn = principle quantum no. n
        pqn = np.asarray(rFblvl['pqn'], 'int64')
        # l is angular moment quantum no. L
        l = rFblvl['l']
        # energy level in inverse cm
        ecm = rFblvl['ecm']
        # statistical weigths/multiplicities
        multr = rFblvl['mult']
        mult = fblvl['mult']
        #
        #
        # for the ionization potential, must use that of the recombined ion
        #
#        iprcm = self.Ipr/const.invCm2Ev
        #
        # get karzas-latter Gaunt factors
        if hasattr(self, 'Klgfb'):
            klgfb = self.Klgfb
        else:
            self.Klgfb = io.klgfbRead()
            klgfb = self.Klgfb
        #
        nTemp = temperature.size
        # statistical weigths/multiplicities
        #
        #
        #wecm=1.e+8/(ipcm-ecm)
        #
        # sometime the rFblvl file does not exist
        if 'mult' in fblvl.keys() and 'mult' in rFblvl.keys():
            fbrate = np.zeros((nlvls,nTemp),np.float64)
            ratg = np.zeros((nlvls),np.float64)
            for ilvl in range(nlvls):
                # scaled energy is relative to the ionization potential of each individual level
                # will add the average energy of a free electron to this to get typical photon energy to
                # evaluate the gaunt factor
                hnuEv = 1.5*const.boltzmann*temperature/const.ev2Erg
                iprLvlEv = self.Ipr - const.invCm2Ev*ecm[ilvl]
                scaledE = np.log(hnuEv/iprLvlEv)
                thisGf = klgfb['klgfb'][pqn[ilvl]-1, l[ilvl]]
                spl = splrep(klgfb['pe'], thisGf)
                gf = np.exp(splev(scaledE, spl))
                ratg[ilvl] = float(multr[ilvl])/float(mult[0]) # ratio of statistical weights
                iprLvlErg = const.ev2Erg*iprLvlEv
                fbrate[ilvl] = ratg[ilvl]*(iprLvlErg**2/float(pqn[ilvl]))*gf/np.sqrt(temperature)
            fbRate = abund*gIoneq*const.freeBoundLoss*(fbrate.sum(axis=0))
        else:
            fbRate = np.zeros((nTemp),np.float64)
        self.FreeBoundLoss = {'rate':fbRate, 'temperature':temperature}

    def freeBoundwB(self, wavelength, includeAbundance=True, includeIoneq=True, useVerner=True, **kwargs):
        """
        Calculate the free-bound emission of an ion. The result is returned as a 2D array to the
        `free_bound_emission` attribute.

        The total free-bound continuum emissivity is given by,

        .. math::
           \\frac{dW}{dtdVd\lambda} = \\frac{1}{4\pi}\\frac{2}{hk_Bc^3m_e\sqrt{2\pi k_Bm_e}}\\frac{E^5}{T^{3/2}}\sum_i\\frac{\omega_i}{\omega_0}\sigma_i^{bf}\exp\left(-\\frac{E - I_i}{k_BT}\\right)

        where :math:`E=hc/\lambda` is the photon energy, :math:`\omega_i` and :math:`\omega_0`
        are the statistical weights of the :math:`i^{\mathrm{th}}` level of the recombined ion
        and the ground level of the recombing ion, respectively, :math:`\sigma_i^{bf}` is the
        photoionization cross-section, and :math:`I_i` is the ionization potential of level :math:`i`.
        This expression comes from Eq. 12 of [3]_. For more information about the free-bound continuum
        calculation, see `Peter Young's notes on free-bound continuum`_.

        The photoionization cross-sections are calculated using the methods of [2]_ for the
        transitions to the ground state and [1]_ for all other transitions. See
        `verner_cross_section` and `karzas_cross_section` for more details.

        .. _Peter Young's notes on free-bound continuum: http://www.pyoung.org/chianti/freebound.pdf

        The free-bound emission is in units of erg
        :math:`\mathrm{cm}^3\mathrm{s}^{-1}\mathrm{\mathring{A}}^{-1}\mathrm{str}^{-1}`. If the emission
        measure has been set, the units will be multiplied by :math:`\mathrm{cm}^{-5}` or
        :math:`\mathrm{cm}^{-3}`, depending on whether it is the line-of-sight or volumetric
        emission measure, respectively.

        Parameters
        ----------
        wavelength : array-like
            In units of angstroms
        include_abundance : `bool`, optional
            If True, include the ion abundance in the final output.
        include_ioneq : `bool`, optional
            If True, include the ionization equilibrium in the final output
        use_verner : `bool`, optional
            If True, cross-sections of ground-state transitions using [2]_, i.e. `verner_cross_section`

        Raises
        ------
        ValueError
            If no .fblvl file is available for this ion

        References
        ----------
        .. [1] Karzas and Latter, 1961, ApJSS, `6, 167
            <http://adsabs.harvard.edu/abs/1961ApJS....6..167K>`_
        .. [2] Verner & Yakovlev, 1995, A&AS, `109, 125
            <http://adsabs.harvard.edu/abs/1995A%26AS..109..125V>`_
        .. [3] Young et al., 2003, ApJSS, `144, 135
            <http://adsabs.harvard.edu/abs/2003ApJS..144..135Y>`_
        """
        wavelength = np.atleast_1d(wavelength)
        if wavelength.size < 2:
            print(' wavelength must have at least two values, current length %3i'%(wavelength.size))
            return
        self.NWavelength = wavelength.size
        # calculate the photon energy in erg
        photon_energy = const.planck*(1.e8*const.light)/wavelength
        prefactor = (2./np.sqrt(2.*np.pi)/(4.*np.pi)/(const.planck*(const.light**3)
                     * (const.emass*const.boltzmann)**(3./2.)))
        # read the free-bound level information for the recombined and recombining ion
        recombining_fblvl = io.fblvlRead(self.IonStr)
        # get the multiplicity of the ground state of the recombining ion
        if 'errorMessage' in recombining_fblvl:
            omega_0 = 1.
        else:
            omega_0 = recombining_fblvl['mult'][0]

        self.Recombined_fblvl = io.fblvlRead(self.nameDict['lower'])
        if 'errorMessage' in self.Recombined_fblvl:
#            raise ValueError('No free-bound information available for {}'.format(util.zion2name(self.Z, self.Stage)))
            errorMessage = 'No free-bound information available for {}'.format(util.zion2name(self.Z, self.Stage))
            fb_emiss = np.zeros((self.NTemperature, self.NWavelength), np.float64)
#            self.free_bound_emission = fb_emiss.squeeze()
            self.FreeBound = {'intensity':fb_emiss, 'temperature':self.Temperature,'wvl':wavelength,'em':self.Em, 'errorMessage':errorMessage}
            return

        energy_over_temp_factor = np.outer(1./(self.Temperature**1.5), photon_energy**5.).squeeze()
#        if self.NWavelength > 1:
#            print(' energy shape %5i %5i'%(energy_over_temp_factor.shape[0],energy_over_temp_factor.shape[1]))
#        else:
#            print(' energy size %5i'%(energy_over_temp_factor.size))
        # sum over levels of the recombined ion
        sum_factor = np.zeros((len(self.Temperature), len(wavelength)))
        for i,omega_i in enumerate(self.Recombined_fblvl['mult']):
            # ionization potential for level i
#            ip = self.ionization_potential - recombined_fblvl['ecm'][i]*const.planck*const.light
            ip = self.IprErg - self.Recombined_fblvl['ecm'][i]*const.planck*const.light
            # skip level if photon energy is not sufficiently high
            if ip < 0. or np.all(np.max(photon_energy) < (self.ionization_potential - ip)):
                continue
            # calculate cross-section
            if i == 0 and useVerner:
                cross_section = self.vernerCross(photon_energy)
            else:
                cross_section = self.karzasCross(photon_energy, ip,
                                                          self.Recombined_fblvl['pqn'][i],
                                                          self.Recombined_fblvl['l'][i])
            scaled_energy = np.outer(1./(const.boltzmann*self.Temperature), photon_energy - ip)
            # the exponential term can go to infinity for low temperatures
            # but if the cross-section is zero this does not matter
            scaled_energy[:,np.where(cross_section == 0.0)] = 0.0
            sum_factor += omega_i/omega_0*np.exp(-scaled_energy)*cross_section

        # combine factors
        fb_emiss = prefactor*energy_over_temp_factor*sum_factor.squeeze()
#        if self.NWavelength > 1:
#            print(' fb emiss.shape %5i %5i'%(fb_emiss.shape[0], fb_emiss.shape[1]))
#        else:
#            print(' fb emiss.size %5i'%(fb_emiss.size))
        # include abundance, ionization equilibrium, photon conversion, emission measure
        if includeAbundance:
            fb_emiss *= self.Abundance
            includeAbundance = self.Abundance
        if includeIoneq:
            if self.NTemperature > 1:
                if self.NWavelength > 1:
#                    fb_emiss *= self.ioneq_one(self.Stage, **kwargs)[:,np.newaxis]
                    fb_emiss *= self.IoneqOne[:,np.newaxis]
                    includeAbundance = self.IoneqOne[:,np.newaxis]
                else:
                    fb_emiss *= self.IoneqOne
                    includeAbundance = self.IoneqOne
            else:
                fb_emiss *= self.IoneqOne
                includeAbundance = self.IoneqOne
        if self.Em is not None:
            if self.Em.size > 1:
                fb_emiss *= self.Em[:,np.newaxis]
            else:
                fb_emiss *= self.Em

        if chdata.Defaults['flux'] == 'photon':
            fb_emiss /= photon_energy
        # the final units should be per angstrom
        fb_emiss /= 1e8

#        self.free_bound_emission = fb_emiss.squeeze()
        self.FreeBound = {'intensity':fb_emiss.squeeze(), 'temperature':self.Temperature,'wvl':wavelength,'em':self.Em, 'ions':self.IonStr,  'abundance':includeAbundance, 'ioneq':includeIoneq}

    def freeBound(self, wvl, verner=1):
        '''
        to calculate the free-bound (radiative recombination) continuum rate coefficient of an ion, where
        the ion is taken to be the target ion,
        including the elemental abundance and the ionization equilibrium population
        uses the Gaunt factors of Karzas, W.J, Latter, R, 1961, ApJS, 6, 167
        for recombination to the ground level, the photoionization cross sections of
        Verner and Yakovlev, 1995, A&ASS, 109, 125
        are used to develop the free-bound cross section
        includes the elemental abundance and the ionization fraction
        provides emissivity = ergs cm^-2 s^-1 str^-1 Angstrom ^-1
        '''
        wvl = np.asarray(wvl, np.float64)
        temperature = self.Temperature
        hnu = 1.e+8*const.planck*const.light/wvl
        #
        if hasattr(self, 'IoneqOne'):
            gIoneq = self.IoneqOne
        else:
            self.ioneqOne()
            gIoneq = self.IoneqOne
        #
        # put in freefree to go through ipymspectrum
        if not np.any(gIoneq) > 0:
            self.FreeBound = {'errorMessage':' no non-zero values of ioneq'}
            return
        #
        em = self.Em
        #
        # the target ion contains that data for fblvl
        #
        if hasattr(self,'Fblvl'):
            fblvl = self.Fblvl
            if 'errorMessage' in fblvl.keys():
                self.FreeBound = fblvl
                return
        elif self.Z == self.Stage-1:
            #dealing with the fully ionized stage
            self.Fblvl = {'mult':[2., 2.]}
            fblvl = self.Fblvl
        else:
            fblvlname = self.nameDict['filename']+'.fblvl'
            if os.path.isfile(fblvlname):
                self.Fblvl = io.fblvlRead(self.IonStr)
                fblvl = self.Fblvl
            # in case there is no fblvl file
            else:
                self.FreeBound = {'errorMessage':' no fblvl file for ion %s'%(self.IonStr)}
                return
        #
        #  need data for the recombined ion
        #
        if hasattr(self,'rFblvl'):
            rfblvl = self.rFblvl
        else:
            lower = self.nameDict['lower']
            lowerDict = util.convertName(lower)
            fblvlname = lowerDict['filename'] +'.fblvl'
            if os.path.isfile(fblvlname):
                self.rFblvl = io.fblvlRead(lower)
                rfblvl = self.rFblvl
            else:
                self.FreeBound = {'errorMessage':' no fblvl file for ion %s'%(self.IonStr)}
                return
        #
        #
        abund = self.Abundance
        #
        #
        nlvls = len(rfblvl['lvl'])
        # pqn = principle quantum no. n
        pqn = rfblvl['pqn']
        # l is angular moment quantum no. L
        l = rfblvl['l']
        # energy level in inverse cm
        ecm = rfblvl['ecm']
        # statistical weigths/multiplicities
        multr = rfblvl['mult']
        mult = fblvl['mult']
        #
        #
        # for the ionization potential, must use that of the recombined ion
        #
        iprcm = self.Ipr/const.invCm2Ev
        #
        # get karzas-latter Gaunt factors
        if hasattr(self,'Klgfb'):
            klgfb = self.Klgfb
        else:
            self.Klgfb = io.klgfbRead()
            klgfb = self.Klgfb
        #
        nWvl = wvl.size
        nTemp = temperature.size
        #
        if verner:
            lvl1 = 1
        else:
            lvl1 = 0
            #
        nWvl = wvl.size
        nTemp = temperature.size
        #
        if verner:
            self.vernerCross(wvl)
            vCross = self.VernerCross
        #
        if (nTemp > 1) and (nWvl > 1):
            mask = np.zeros((nlvls,nTemp,nWvl),np.bool_)
            fbrate = np.zeros((nlvls,nTemp,nWvl),np.float64)
            fbRate = np.zeros((nTemp,nWvl),np.float64)
            expf = np.zeros((nlvls,nTemp,nWvl),np.float64)
            ratg = np.zeros((nlvls),np.float64)
            ratg[0] = float(multr[0])/float(mult[0])
            iprLvlEv = self.Ipr - const.invCm2Ev*ecm[0]
            iprLvlErg = const.ev2Erg*iprLvlEv
            iprLvlCm = (iprcm - ecm[0])
            for itemp in range(nTemp):
                mask[0,itemp] = 1.e+8/wvl < (iprcm - ecm[0])
                expf[0,itemp] = np.exp((iprLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature[itemp]))
                fbrate[0,itemp] = em[itemp]*abund*gIoneq[itemp]*(const.planck*const.light/(1.e-8*wvl))**5*const.verner*ratg[0]*expf[0,itemp]*vCross/temperature[itemp]**1.5
            for ilvl in range(lvl1,nlvls):
                iprLvlEv = self.Ipr - const.invCm2Ev*ecm[ilvl]
                iprLvlErg = const.ev2Erg*iprLvlEv
                scaledE = np.log(const.ev2Ang/(iprLvlEv*wvl))
                thisGf = klgfb['klgfb'][pqn[ilvl]-1, l[ilvl]]
                spl = splrep(klgfb['pe'], thisGf)
                gf = np.exp(splev(scaledE, spl))
                ratg[ilvl] = float(multr[ilvl])/float(mult[0]) # ratio of statistical weights
            #
                for itemp in range(nTemp):
                    expf[ilvl] = np.exp((iprLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature[itemp]))
                    expf[ilvl,itemp] = np.exp((iprLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature[itemp]))
                    mask[ilvl,itemp] = 1.e+8/wvl < (iprcm - ecm[ilvl])
                    fbrate[ilvl,itemp] = em[itemp]*abund*gIoneq[itemp]*const.freeBound*ratg[ilvl]*(iprLvlErg**2/float(pqn[ilvl]))*gf*expf[ilvl,itemp]/(temperature[itemp]**1.5*(wvl)**2)
            fbrma = np.ma.array(fbrate)
            fbrma.mask =  mask
            fbrma.fill_value = 0.
            fbIntensity = fbrma.sum(axis=0)
#            for itemp in range(nTemp):
#                fbRate += em[itemp]*abund*gIoneq[itemp]*fbrma[itemp]
#            fbRate = fbrma.sum(axis=0)
#            fbRate.fill_value = 0.
            self.FreeBound = {'intensity':fbIntensity, 'temperature':temperature,'wvl':wvl,'em':em}
            #
        elif (nTemp == 1) and (nWvl > 1):
            mask = np.zeros((nlvls,nWvl),np.bool_)
            fbrate = np.zeros((nlvls,nWvl),np.float64)
            expf = np.zeros((nlvls,nWvl),np.float64)
            ratg = np.zeros((nlvls),np.float64)
            # mask is true for bad values
            ratg[0] = float(multr[0])/float(mult[0])
            iprLvlEv = self.Ipr - const.invCm2Ev*ecm[0]
            iprLvlErg = const.ev2Erg*iprLvlEv
            iprLvlCm = (iprcm - ecm[0])
            #
            mask[0] = 1.e+8/wvl < iprcm
            expf[0] = np.exp((iprLvlErg - hnu)/(const.boltzmann*temperature))
            # both expressions for fbrate[0] match the IDL output
            fbrate[0] = (const.planck*const.light/(1.e-8*wvl))**5*const.verner*ratg[0]*expf[0]*vCross/temperature**1.5
            # factor of 1.e-8 converts to Angstrom^-1, otherwise it would be cm^-1
#            fbrate[0] = 1.e-8*const.freeBounde*hnu**5*ratg[0]*expf[0]*vCross/temperature**1.5
            #
            for ilvl in range(lvl1,nlvls):
                iprLvlEv = self.Ipr - const.invCm2Ev*ecm[ilvl]
                iprLvlErg = const.ev2Erg*iprLvlEv
                iprLvlCm = (iprcm - ecm[ilvl])
                # scaled energy is relative to the ionization potential of each individual level
                scaledE = np.log(const.ev2Ang/(iprLvlEv*wvl))
                thisGf = klgfb['klgfb'][pqn[ilvl]-1, l[ilvl]]
                spl = splrep(klgfb['pe'], thisGf)
                gf = np.exp(splev(scaledE, spl))
                mask[ilvl] = 1.e+8/wvl < iprLvlCm
                ratg[ilvl] = float(multr[ilvl])/float(mult[0]) # ratio of statistical weights
                expf[ilvl] = np.exp((iprLvlErg - hnu)/(const.boltzmann*temperature))
                fbrate[ilvl] = const.freeBound*ratg[ilvl]*(iprLvlErg**2/float(pqn[ilvl]))*expf[ilvl]*gf/(temperature**1.5*(wvl)**2)
            fbrma = np.ma.array(fbrate)
            fbrma.mask =  mask
            fbrma.fill_value = 0.
            fbRate = em*abund*gIoneq*fbrma.sum(axis=0)
            fbRate.fill_value = 0.
            self.FreeBound = {'fbRate':fbRate, 'intensity':fbRate.data, 'temperature':temperature,'wvl':wvl, 'mask':mask, 'expf':expf,'vCross':vCross}
        #elif (nTemp > 1) and (nWvl == 1):
        else:
            self.FreeBound = {'intensity':np.zeros(nTemp,np.float64),'errorMessage':' this is the case of a single wavelength'}

    def freeBoundEmiss(self, wvl, verner=1):

        """
        Calculates the free-bound (radiative recombination) continuum emissivity of an ion.
        Provides emissivity in units of ergs :math:`\mathrm{cm}^{-2}` :math:`\mathrm{s}^{-1}` :math:`\mathrm{str}^{-1}` :math:`\mathrm{\AA}^{-1}` for an individual ion.

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
        wvl = np.asarray(wvl, np.float64)
        temperature = self.Temperature
        hnu = 1.e+8*const.planck*const.light/wvl
        #
        #
        em = self.Em
        #
        # the target ion contains that data for fblvl
        #
        if hasattr(self,'Fblvl'):
            fblvl = self.Fblvl
            if 'errorMessage' in fblvl.keys():
                self.FreeBound = fblvl
                return
        elif self.Z == self.Stage-1:
            #dealing with the fully ionized stage
            self.Fblvl = {'mult':[2., 2.]}
            fblvl = self.Fblvl
        else:
            fblvlname = self.nameDict['filename']+'.fblvl'
            if os.path.isfile(fblvlname):
                self.Fblvl = io.fblvlRead(self.IonStr)
                fblvl = self.Fblvl
            # in case there is no fblvl file
            else:
                self.FreeBound = {'errorMessage':' no fblvl file for ion %s'%(self.IonStr)}
                return
        #
        #  need data for the recombined ion
        #
        if hasattr(self,'rFblvl'):
            rfblvl = self.rFblvl
        else:
            lower = self.nameDict['lower']
            lowerDict = util.convertName(lower)
            fblvlname = lowerDict['filename'] +'.fblvl'
            if os.path.isfile(fblvlname):
                self.rFblvl = io.fblvlRead(lower)
                rfblvl = self.rFblvl
            else:
                self.FreeBound = {'errorMessage':' no fblvl file for ion %s'%(self.IonStr)}
                return
        #
        #
        nlvls = len(rfblvl['lvl'])
        # pqn = principle quantum no. n
        pqn = rfblvl['pqn']
        # l is angular moment quantum no. L
        l = rfblvl['l']
        # energy level in inverse cm
        ecm = rfblvl['ecm']
        # statistical weigths/multiplicities
        multr = rfblvl['mult']
        mult = fblvl['mult']
        #
        #
        # for the ionization potential, must use that of the recombined ion
        #
        iprcm = self.Ipr/const.invCm2Ev
        #
        # get karzas-latter Gaunt factors
        if hasattr(self,'Klgfb'):
            klgfb = self.Klgfb
        else:
            self.Klgfb = io.klgfbRead()
            klgfb = self.Klgfb
        #
        nWvl = wvl.size
        nTemp = temperature.size
        #
        if verner:
            lvl1 = 1
        else:
            lvl1 = 0
            #
        nWvl = wvl.size
        nTemp = temperature.size
        #
        if verner:
            self.vernerCross(wvl)
            vCross = self.VernerCross
        #
        if (nTemp > 1) and (nWvl > 1):
            mask = np.zeros((nlvls,nTemp,nWvl),np.bool_)
            fbrate = np.zeros((nlvls,nTemp,nWvl),np.float64)
            fbRate = np.zeros((nTemp,nWvl),np.float64)
            expf = np.zeros((nlvls,nTemp,nWvl),np.float64)
            ratg = np.zeros((nlvls),np.float64)
            ratg[0] = float(multr[0])/float(mult[0])
            iprLvlEv = self.Ipr - const.invCm2Ev*ecm[0]
            iprLvlErg = const.ev2Erg*iprLvlEv
            iprLvlCm = (iprcm - ecm[0])
            for itemp in range(nTemp):
                mask[0,itemp] = 1.e+8/wvl < (iprcm - ecm[0])
                expf[0,itemp] = np.exp((iprLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature[itemp]))
                fbrate[0,itemp] = em[itemp]*(const.planck*const.light/(1.e-8*wvl))**5*const.verner*ratg[0]*expf[0,itemp]*vCross/temperature[itemp]**1.5
            for ilvl in range(lvl1,nlvls):
                iprLvlEv = self.Ipr - const.invCm2Ev*ecm[ilvl]
                iprLvlErg = const.ev2Erg*iprLvlEv
                scaledE = np.log(const.ev2Ang/(iprLvlEv*wvl))
                thisGf = klgfb['klgfb'][pqn[ilvl]-1, l[ilvl]]
                spl = splrep(klgfb['pe'], thisGf)
                gf = np.exp(splev(scaledE, spl))
                ratg[ilvl] = float(multr[ilvl])/float(mult[0]) # ratio of statistical weights
            #
                for itemp in range(nTemp):
                    expf[ilvl] = np.exp((iprLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature[itemp]))
                    expf[ilvl,itemp] = np.exp((iprLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature[itemp]))
                    mask[ilvl,itemp] = 1.e+8/wvl < (iprcm - ecm[ilvl])
                    fbrate[ilvl,itemp] = em[itemp]*const.freeBound*ratg[ilvl]*(iprLvlErg**2/float(pqn[ilvl]))*gf*expf[ilvl,itemp]/(temperature[itemp]**1.5*(wvl)**2)
            fbrma = np.ma.array(fbrate)
            fbrma.mask =  mask
            fbrma.fill_value = 0.
            fbIntensity = fbrma.sum(axis=0)
#            for itemp in range(nTemp):
#                fbRate += em[itemp]*abund*gIoneq[itemp]*fbrma[itemp]
#            fbRate = fbrma.sum(axis=0)
#            fbRate.fill_value = 0.
            self.FreeBoundEmiss = {'emiss':fbIntensity, 'temperature':temperature,'wvl':wvl,'em':em}
            #
        elif (nTemp == 1) and (nWvl > 1):
            mask = np.zeros((nlvls,nWvl),np.bool_)
            fbrate = np.zeros((nlvls,nWvl),np.float64)
            expf = np.zeros((nlvls,nWvl),np.float64)
            ratg = np.zeros((nlvls),np.float64)
            # mask is true for bad values
            ratg[0] = float(multr[0])/float(mult[0])
            iprLvlEv = self.Ipr - const.invCm2Ev*ecm[0]
            iprLvlErg = const.ev2Erg*iprLvlEv
            iprLvlCm = (iprcm - ecm[0])
            #
            mask[0] = 1.e+8/wvl < iprcm
            expf[0] = np.exp((iprLvlErg - hnu)/(const.boltzmann*temperature))
            # both expressions for fbrate[0] match the IDL output
            fbrate[0] = (const.planck*const.light/(1.e-8*wvl))**5*const.verner*ratg[0]*expf[0]*vCross/temperature**1.5
            # factor of 1.e-8 converts to Angstrom^-1, otherwise it would be cm^-1
#            fbrate[0] = 1.e-8*const.freeBounde*hnu**5*ratg[0]*expf[0]*vCross/temperature**1.5
            #
            for ilvl in range(lvl1,nlvls):
                iprLvlEv = self.Ipr - const.invCm2Ev*ecm[ilvl]
                iprLvlErg = const.ev2Erg*iprLvlEv
                iprLvlCm = (iprcm - ecm[ilvl])
                # scaled energy is relative to the ionization potential of each individual level
                scaledE = np.log(const.ev2Ang/(iprLvlEv*wvl))
                thisGf = klgfb['klgfb'][pqn[ilvl]-1, l[ilvl]]
                spl = splrep(klgfb['pe'], thisGf)
                gf = np.exp(splev(scaledE, spl))
                mask[ilvl] = 1.e+8/wvl < iprLvlCm
                ratg[ilvl] = float(multr[ilvl])/float(mult[0]) # ratio of statistical weights
                expf[ilvl] = np.exp((iprLvlErg - hnu)/(const.boltzmann*temperature))
                fbrate[ilvl] = const.freeBound*ratg[ilvl]*(iprLvlErg**2/float(pqn[ilvl]))*expf[ilvl]*gf/(temperature**1.5*(wvl)**2)
            fbrma = np.ma.array(fbrate)
            fbrma.mask =  mask
            fbrma.fill_value = 0.
            fbRate = em*fbrma.sum(axis=0)
            fbRate.fill_value = 0.
            self.FreeBoundEmiss = {'emiss':fbRate.data, 'temperature':temperature,'wvl':wvl, 'em':em}
        #elif (nTemp > 1) and (nWvl == 1):
        else:
            self.FreeBoundEmiss = {'emiss':np.zeros(nTemp,np.float64),'errorMessage':' this is the case of a single wavelength'}


    def vernerCross(self, wvl):
        """
        Calculates the photoionization cross-section using data from [1]_ for
        transitions to the ground state.

        The photoionization cross-section can be expressed as :math:`\sigma_i^{fb}=F(E/E_0)` where
        :math:`F` is an analytic fitting formula given by Eq. 1 of [1]_,

        .. math::
           F(y) = ((y-1)^2 + y_w^2)y^{-Q}(1 + \sqrt{y/y_a})^{-P},

        where :math:`E` is the photon energy, :math:`n` is the principal quantum number,
        :math:`l` is the orbital quantum number, :math:`Q = 5.5 + l - 0.5P`, and
        :math:`\sigma_0,E_0,y_w,y_a,P` are fitting paramters. These can be read in using
        `ChiantiPy.tools.io.vernerRead`.

        References
        ----------
        .. [1] Verner & Yakovlev, 1995, A&AS, `109, 125
            <http://adsabs.harvard.edu/abs/1995A%26AS..109..125V>`_
        """
        # read verner data
        verner_info = io.vernerRead()
        eth = verner_info['eth'][self.Z,self.Stage-1]   #*const.ev2Erg
        yw = verner_info['yw'][self.Z,self.Stage-1]
        ya = verner_info['ya'][self.Z,self.Stage-1]
        p = verner_info['p'][self.Z,self.Stage-1]

        # convert from megabarn to cm^2
        sigma0 = verner_info['sig0'][self.Z,self.Stage-1]*1e-18
        e0 = verner_info['e0'][self.Z,self.Stage-1]  #*const.ev2Erg
        q = 5.5 + verner_info['l'][self.Z,self.Stage-1] - 0.5*p

        # scaled photon energy
        en = const.ev2Ang/wvl
        y = en/e0
        # fitting function
        F = ((y - 1.)**2 + yw**2)*(y**(-q))*(1. + np.sqrt(y/ya))**(-p)
        cross_section = sigma0*F

        self.VernerCross = np.where(en < eth, 0., cross_section)

    def karzasCross(self, photon_energy, ionization_potential, n, l):
        """
        Calculate the photoionization cross-sections using the Gaunt factors of [1]_.

        The free-bound photoionization cross-section is given by,

        .. math::
           \sigma_i^{bf} = 1.077294\\times8065.54\\times10^{16}\left(\\frac{I_i}{hc}\\right)^2\left(\\frac{hc}{E}\\right)^3\\frac{g_{bf}}{n_i},

        where :math:`I_i` is the ionization potential of the :math:`i^{\mathrm{th}}` level,
        :math:`E` is the photon energy, :math:`g_{bf}` is the Gaunt factor calculated
        according to [1]_, and :math:`n_i` is the principal quantum number of the
        :math:`i^{\mathrm{th}}` level. :math:`\sigma_i^{bf}` is units of :math:`\mathrm{cm}^{2}`.
        This expression is given by Eq. 13 of [2]_. For more information on the photoionization
        cross-sections, see `Peter Young's notes on free-bound continuum`_.

        .. _Peter Young's notes on free-bound continuum: http://www.pyoung.org/chianti/freebound.pdf

        Parameters
        ----------
        photon_energy : array-like
        ionization_potential : `float`
        n : `int`
        l : `int`

        References
        ----------
        .. [1] Karzas and Latter, 1961, ApJSS, `6, 167
            <http://adsabs.harvard.edu/abs/1961ApJS....6..167K>`_
        .. [2] Young et al., 2003, ApJSS, `144, 135
            <http://adsabs.harvard.edu/abs/2003ApJS..144..135Y>`_
        """
        # numerical constant, in Mbarn
        kl_constant = 1.077294e-1*8065.54e3
        # read in KL gaunt factor data
        karzas_info = io.klgfbRead()
        if n <= karzas_info['klgfb'].shape[0]:
            scaled_energy = np.log10(photon_energy/ionization_potential)
            f_gf = splrep(karzas_info['pe'], karzas_info['klgfb'][n-1,l,:])
            gaunt_factor = np.exp(splev(scaled_energy, f_gf))
        else:
            gaunt_factor = 1.

        # scaled energy factor, converted to cm^-1
        energy_factor = (((ionization_potential/const.planck/const.light)**2.)
                         * ((photon_energy/const.planck/const.light)**(-3)))
        # cross-section, convert to cm^2
        cross_section = (kl_constant*energy_factor*gaunt_factor/n)*1e-18

        return np.where(photon_energy >= ionization_potential, cross_section, 0.)



    def klgfbInterp(self, wvl, n, l):
        '''A Python version of the CHIANTI IDL procedure karzas_xs.

        Interpolates free-bound gaunt factor of Karzas and Latter, (1961, Astrophysical Journal
        Supplement Series, 6, 167) as a function of wavelength (wvl).'''
        try:
            klgfb = self.Klgfb
        except:
            self.Klgfb = util.klgfbRead()
            klgfb = self.Klgfb
        # get log of photon energy relative to the ionization potential
        sclE = np.log(self.Ip/(wvl*const.ev2ang))
        thisGf = klgfb['klgfb'][n-1, l]
        spl = splrep(klgfb['pe'], thisGf)
        gf = splev(sclE, spl)
        return gf

    def ioneqOne(self):
        '''
        Provide the ionization equilibrium for the selected ion as a function of temperature.
        Similar to but not identical to ion.ioneqOne() - the ion class needs to be able to handle
        the 'dielectronic' ions
        returned in self.IoneqOne
        '''
        #
        if hasattr(self, 'Temperature'):
            temperature = self.Temperature
        else:
            return
        #
        if hasattr(self, 'IoneqAll'):
            ioneqAll = self.IoneqAll
        else:
            self.IoneqAll = io.ioneqRead(ioneqName = self.Defaults['ioneqfile'])
            ioneqAll = self.IoneqAll
        #
        ioneqTemperature = ioneqAll['ioneqTemperature']
        Z = self.Z
        stage = self.Stage
        ioneqOne = np.zeros_like(temperature)
        #
        thisIoneq = ioneqAll['ioneqAll'][Z-1,stage-1].squeeze()
        gioneq = thisIoneq > 0.
        goodt1 = self.Temperature >= ioneqTemperature[gioneq].min()
        goodt2 = self.Temperature <= ioneqTemperature[gioneq].max()
        goodt = np.logical_and(goodt1,goodt2)
        y2 = splrep(np.log(ioneqTemperature[gioneq]),np.log(thisIoneq[gioneq]),s=0)
        #
        if goodt.sum() > 0:
            if self.Temperature.size > 1:
                gIoneq = splev(np.log(self.Temperature[goodt]),y2)   #,der=0)
                ioneqOne[goodt] = np.exp(gIoneq)
            else:
                gIoneq = splev(np.log(self.Temperature),y2)
                ioneqOne = np.exp(gIoneq)
                ioneqOne = np.atleast_1d(ioneqOne)
            self.IoneqOne = ioneqOne
        else:
            self.IoneqOne = np.zeros_like(self.Temperature)


    def ioneq_one(self, stage, **kwargs):
        """
        Calculate the equilibrium fractional ionization of the ion as a function of temperature.

        Uses the `ChiantiPy.core.ioneq` module and does a first-order spline interpolation to the data. An
        ionization equilibrium file can be passed as a keyword argument, `ioneqfile`. This can
        be passed through as a keyword argument to any of the functions that uses the
        ionization equilibrium.

        Parameters
        ----------
        stage : int
            Ionization stage, e.g. 25 for Fe XXV
        """
        tmp = ioneq(self.Z)
        tmp.load(ioneqName=kwargs.get('ioneqfile', None))
        ionization_equilibrium = splev(self.Temperature,
                                       splrep(tmp.Temperature, tmp.Ioneq[stage-1,:], k=1), ext=1)
        return np.where(ionization_equilibrium < 0., 0., ionization_equilibrium)
