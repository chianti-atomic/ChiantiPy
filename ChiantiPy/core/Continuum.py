"""
Continuum module

"""

import os

import numpy as np
from scipy.interpolate import splev, splrep
from scipy import special
from scipy.ndimage import map_coordinates
#import matplotlib.pyplot as plt
#import ChiantiPy.core as ch
from .Ioneq import ioneq
from ChiantiPy.base import ionTrails
import ChiantiPy.tools.data as chdata
import ChiantiPy.tools.util as util
import ChiantiPy.tools.io as io
import ChiantiPy.tools.constants as const


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

    Keyword arguments
    -----------------
    abundance : `float`, `str`, optional
        Elemental abundance relative to Hydrogen or name of CHIANTI abundance file,
        without the '.abund' suffix, e.g. 'sun_photospheric_1998_grevesse'.
    em : array-like, optional
        Line-of-sight emission measure (:math:`\int\mathrm{d}l\,n_en_H`), in units of
        :math:`\mathrm{cm}^{-5}`, or the volumetric emission measure (:math:`\int\mathrm{d}V\,n_en_H`)
        in units of :math:`\mathrm{cm}^{-3}`.
    verbose : `bool`
        if True, prints additional info to the console

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

    References
    ----------

    .. [101] Sutherland, R. S., 1998, MNRAS, `300, 321 <http://adsabs.harvard.edu/abs/1998MNRAS.300..321S>`_
    .. [102] Verner & Yakovlev, 1995, A&AS, `109, 125 <http://adsabs.harvard.edu/abs/1995A%26AS..109..125V>`_
    .. [103] Karzas and Latter, 1961, ApJSS, `6, 167 <http://adsabs.harvard.edu/abs/1961ApJS....6..167K>`_
    .. [104] Itoh, N. et al., 2000, ApJS, `128, 125 <http://adsabs.harvard.edu/abs/2000ApJS..128..125I>`_
    .. [106] Mewe, R. et al., 1986, A&AS, `65, 511 <http://adsabs.harvard.edu/abs/1986A%26AS...65..511M>`_
    .. [107] Rybicki and Lightman, 1979, Radiative Processes in Astrophysics,
            `(Wiley-VCH) <http://adsabs.harvard.edu/abs/1986rpa..book.....R>`_
    .. [108] Gronenschild, E.H.B.M. and Mewe, R., 1978, A&AS, `32, 283 <http://adsabs.harvard.edu/abs/1978A%26AS...32..283G>`_
    """

    def __init__(self, ionStr,  temperature, abundance=None, em=None, verbose=0):


        self.IonStr = ionStr
        self.nameDict = util.convertName(ionStr)
        #  Z is the ion charge
        self.Z = self.nameDict['Z']
        # Zion is the charge of the ion
        self.Zion = self.nameDict['Ion'] -1
        self.Stage = self.nameDict['Ion']
        self.Ion = self.nameDict['Ion']
        self.Dielectronic = self.nameDict['Dielectronic']
        self.Temperature = np.atleast_1d(temperature)
        self.NTemperature = self.Temperature.size
        self.Defaults = chdata.Defaults

        if self.Defaults['wavelength'] != 'angstrom':
            print(' the continuum can only be calculated for wavelengths in angstroms')
            return

        self.argCheck(temperature=temperature, eDensity=None, pDensity=None, em=em, verbose=verbose)

        self.Labels = util.units(chdata.Defaults)


        self.Ip = chdata.Ip[self.Z-1, self.Stage-1]
        self.Ipr = chdata.Ip[self.Z-1, self.Stage-2]
        self.ionization_potential = chdata.Ip[self.Z-1, self.Stage-1]*const.ev2Erg
        self.IprErg = self.Ipr*const.ev2Erg
        self.IprCm = 1.e+8/(const.ev2Ang/self.Ipr)
        # Set abundance
        if abundance is not None:
            if isinstance(abundance,  float):
                self.Abundance = abundance
                # need the following for p2eRatio but will use default abundance
                self.AbundAll = chdata.AbundanceDefault['abundance']
            elif isinstance(abundance, str):
                ab = io.abundanceRead(abundance,  verbose=1)
                self.Abundance = ab['abundance'][self.Z-1]
                self.AbundanceName = abundance
                self.AbundAll = ab['abundance']

        else:
            self.Abundance = chdata.AbundanceDefault['abundance'][self.Z-1]
            self.AbundAll = chdata.AbundanceDefault['abundance']
        self.ioneqOne()

    def free_free_loss(self,  includeAbund=True, includeIoneq=True, **kwargs):
        """
        Calculate the free-free energy loss rate of an ion. The result is returned to the
        `free_free_loss` attribute.

        The free-free radiative loss rate is given by Eq. 5.15a of [107]_. Writing the numerical
        constant in terms of the fine structure constant :math:`\\alpha`,

        .. math::
           \\frac{dW}{dtdV} = \\frac{4\\alpha^3h^2}{3\pi^2m_e}\left(\\frac{2\pi k_B}{3m_e}\\right)^{1/2}Z^2T^{1/2}\\bar{g}_B

        where where :math:`Z` is the nuclear charge, :math:`T` is the electron temperature, and
        :math:`\\bar{g}_{B}` is the wavelength-averaged and velocity-averaged Gaunt factor. The
        Gaunt factor is calculated using the methods of [103]_. Note that this expression for the
        loss rate is just the integral over wavelength of Eq. 5.14a of [107]_, the free-free emission, and
        is expressed in units of erg :math:`\mathrm{cm}^3\,\mathrm{s}^{-1}`.

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
        if includeAbund:
            prefactor *= self.Abundance
        if includeIoneq:
            prefactor *= self.ioneq_one(self.Stage, **kwargs)

        self.FreeFreeLoss = prefactor*(self.Zion**2)*np.sqrt(self.Temperature)*gaunt_factor

    def freeFree(self, wavelength, includeAbund=True, includeIoneq=True, **kwargs):
        """
        Calculates the free-free emission for a single ion. The result is returned as a dict to
        the `FreeFree` attribute.  The dict has the keywords `intensity`, `wvl`, `temperature`, `em`.

        The free-free emission for the given ion is calculated according Eq. 5.14a of [107]_,
        substituting :math:`\\nu=c/\lambda`, dividing by the solid angle, and writing the numerical
        constant in terms of the fine structure constant :math:`\\alpha`,

        .. math::
           \\frac{dW}{dtdVd\lambda} = \\frac{c}{3m_e}\left(\\frac{\\alpha h}{\pi}\\right)^3\left(\\frac{2\pi}{3m_ek_B}\\right)^{1/2}\\frac{Z^2}{\lambda^2T^{1/2}}\exp{\left(-\\frac{hc}{\lambda k_BT}\\right)}\\bar{g}_{ff},

        where :math:`Z` is the nuclear charge, :math:`T` is the electron temperature in K, and
        :math:`\\bar{g}_{ff}` is the velocity-averaged Gaunt factor. The Gaunt factor is estimated
        using the methods of [104]_ and [101]_, depending on the temperature and energy regime. See
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
        includeAbund : `bool`, optional
            If True, include the ion abundance in the final output.
        includeIoneq : `bool`, optional
            If True, include the ionization equilibrium in the final output

        """
        if self.Defaults['wavelength'] != 'angstrom':
            print(' the continuum can only be calculated for wavelengths in angstroms')
            return


        xlabel = self.Labels['xlabel']
        ylabel = self.Labels['spectrumYlabel']

        if np.array_equal(self.Em, np.ones_like(self.Em)):
            ylabel += '($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'

        wavelength = np.atleast_1d(wavelength)
        # define the numerical prefactor
        prefactor = ((const.light*1e8)/3./const.emass
                     * (const.fine*const.planck/np.pi)**3
                     * np.sqrt(2.*np.pi/3./const.emass/const.boltzmann))
        # include temperature dependence
        prefactor *= self.Zion**2/np.sqrt(self.Temperature)
        if includeAbund:
            prefactor *= self.Abundance
        if includeIoneq:
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
        self.FreeFree = {'intensity':free_free_emission, 'temperature':self.Temperature, 'wvl':wavelength,
            'em':self.Em, 'ions':self.IonStr, 'xlabel':xlabel,  'ylabel':ylabel,  'gf':gf}

    def freeFreeLoss(self, includeAbund=True, includeIoneq=True,  **kwargs):
        """
        Calculate the free-free energy loss rate of an ion. The result is returned to the
        `FreeFreeLoss` attribute.

        The free-free radiative loss rate is given by Eq. 5.15a of [107]_. Writing the numerical
        constant in terms of the fine structure constant :math:`\\alpha`,

        .. math::
           \\frac{dW}{dtdV} = \\frac{4\\alpha^3h^2}{3\pi^2m_e}\left(\\frac{2\pi k_B}{3m_e}\\right)^{1/2}Z^2T^{1/2}\\bar{g}_B

        where where :math:`Z` is the ion charge, :math:`T` is the electron temperature, and
        :math:`\\bar{g}_{B}` is the wavelength-averaged and velocity-averaged Gaunt factor. The
        Gaunt factor is calculated using the methods of [103]_. Note that this expression for the
        loss rate is just the integral over wavelength of Eq. 5.14a of [103]_, the free-free emission, and
        is expressed in units of erg :math:`\mathrm{cm}^3\,\mathrm{s}^{-1}`.

        """
        # interpolate wavelength-averaged  gaunt factors of Sutherland
        gf_kl_info = io.gffintRead()
        gamma_squared = self.IprErg/const.boltzmann/self.Temperature
#        for i, atemp in enumerate(self.Temperature):
#            print('%s T:  %10.2e gamma_squared  %10.2e'%(self.IonStr, atemp, gamma_squared[i]))
        gaunt_factor = splev(np.log(gamma_squared),
                             splrep(gf_kl_info['g2'],gf_kl_info['gffint']), ext=3)
        # calculate numerical constant = Sutherland's Fk
        prefactor = (4.*(const.fine**3)*(const.planck**2)/3./(np.pi**2)/const.emass
                     * np.sqrt(2.*np.pi*const.boltzmann/3./const.emass))
        # include abundance and ionization equilibrium
        if includeAbund:
            prefactor *= self.Abundance
        if includeIoneq:
            prefactor *= self.IoneqOne

        rate = prefactor*(self.Zion**2)*np.sqrt(self.Temperature)*gaunt_factor
        self.FreeFreeLoss = {'temperature':self.Temperature, 'rate':rate,  'gf':gaunt_factor,
            'prefactor':prefactor}


    def itoh_gaunt_factor(self, wavelength):
        """
        Calculates the free-free gaunt factors of [104]_.

        An analytic fitting formulae for the relativistic Gaunt factor is given by Eq. 4 of [104]_,

        .. math::
           g_{Z} = \sum^{10}_{i,j=0}a_{ij}t^iU^j

        where,

        .. math::
           t = \\frac{1}{1.25}(\log_{10}{T} - 7.25),\\
           U = \\frac{1}{2.5}(\log_{10}{u} + 1.5),

        :math:`u=hc/\lambda k_BT`, and :math:`a_{ij}` are the fitting coefficients and are read
        in using `ChiantiPy.tools.io.itohRead` and are given in Table 4 of [104]_. These values
        are valid for :math:`6<\log_{10}(T)< 8.5` and :math:`-4<\log_{10}(u)<1`.

        See Also
        --------
        ChiantiPy.tools.io.itohRead : Read in Gaunt factor coefficients from [104]_

        """
        # calculate scaled energy and temperature
        lower_u = const.planck*(1.e8*const.light)/const.boltzmann/np.outer(self.Temperature, wavelength)
        upper_u = 1./2.5*(np.log10(lower_u) + 1.5)
        t = 1./1.25*(np.log10(self.Temperature) - 7.25)
        # read in Itoh coefficients
        itoh_coefficients = io.itohRead()['itohCoef'][self.Zion - 1].reshape(11,11)
        # calculate Gaunt factor
        gf = np.zeros(upper_u.shape)
        if self.Ntemp == 1:
            for j in range(11):
                for i in range(11):
                    gf += (itoh_coefficients[i,j]*(t**i))*(upper_u**j)
        else:
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
        Calculates the free-free gaunt factor calculations of [101]_.

        The Gaunt factors of [101]_ are read in using `ChiantiPy.tools.io.gffRead`
        as a function of :math:`u` and :math:`\gamma^2`. The data are interpolated
        to the appropriate wavelength and temperature values using
        `~scipy.ndimage.map_coordinates`.

        """
        # calculate scaled quantities
        lower_u = const.planck*(1.e8*const.light)/const.boltzmann/np.outer(self.Temperature,wavelength)
        gamma_squared = (self.Zion**2)*const.ryd2erg/const.boltzmann/self.Temperature[:,np.newaxis]*np.ones(lower_u.shape)
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
        approach of [108]_ and [106]_. Eq. 1a of [106]_ can be integrated over wavelength to get the free-bound loss rate,

        .. math::
           \\frac{dW}{dtdV} = C_{ff}\\frac{k}{hc}T^{1/2}G_{fb},

        in units of erg :math:`\mathrm{cm}^3\,\mathrm{s}^{-1}` where :math:`G_{fb}` is the free-bound Gaunt factor as
        given by Eq. 15 of [106]_ (see `mewe_gaunt_factor` for more details) and :math:`C_{ff}` is the numerical constant
        as given in Eq. 4 of [108]_ and can be written in terms of the fine structure constant :math:`\\alpha`,

        .. math::
           C_{ff}\\frac{k}{hc} = \\frac{8}{3}\left(\\frac{\pi}{6}\\right)^{1/2}\\frac{h^2\\alpha^3}{\pi^2}\\frac{k_B}{m_e^{3/2}} \\approx 1.43\\times10^{-27}

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
        approach of [108]_ and [106]_. Eq. 1a of [106]_ can be integrated over wavelength to get the free-bound loss rate,

        .. math::
           \\frac{dW}{dtdV} = C_{ff}\\frac{k}{hc}T^{1/2}G_{fb},

        in units of erg :math:`\mathrm{cm}^3\,\mathrm{s}^{-1}` where :math:`G_{fb}` is the free-bound Gaunt factor as
        given by Eq. 15 of [106]_ (see `mewe_gaunt_factor` for more details) and :math:`C_{ff}` is the numerical constant
        as given in Eq. 4 of [108]_ and can be written in terms of the fine structure constant :math:`\\alpha`,

        .. math::
           C_{ff}\\frac{k}{hc} = \\frac{8}{3}\left(\\frac{\pi}{6}\\right)^{1/2}\\frac{h^2\\alpha^3}{\pi^2}\\frac{k_B}{m_e^{3/2}} \\approx 1.43\\times10^{-27}

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

    def freeBoundLossMao(self,  includeAbund=False,  includeIoneq=False):
        """
        to calculate the radiative loss rate from the parameters of
        Mao J., Kaastra J., Badnell N.R.
        <Astron. Astrophys. 599, A10 (2017)>=2017A&A...599A..10M

        """
        pars = io.maoParsRead()
        nameDict = util.convertName(self.IonStr)
        Z = nameDict['Z']
        Iso = nameDict['iso']
        ionIndex = 1000*(Iso+1) + Z
#        print(' %s  Z %i  Iso %i'%(self.IonStr, Z,  Iso))
        Index = []
        for iso,  z in zip(pars['seq'], pars['z']):
            Index.append(1000*(iso) + z)
        idx = Index.index(ionIndex)
        #
        abund = self.Abundance
        ioneq = self.IoneqOne
        em = self.Em
        temp = self.Temperature
        tev = const.boltzmannEv*self.Temperature
        a0 = pars['a0'][idx]
        a1 = pars['a1'][idx]
        a2 = pars['a2'][idx]
        b0 = pars['b0'][idx]
        b1 = pars['b1'][idx]
        b2 = pars['b2'][idx]
        c0 = pars['c0'][idx]
#        print('a0:  %10.2e b0:  %10.2e  c0:  %10.2e'%(a0,  b0, c0))
#        print('a0:  %10.2e a1:  %10.2e  a2:  %10.2e'%(a0,  a1, a2))
#        print('b0:  %10.2e b1:  %10.2e  b2:  %10.2e'%(b0,  b1, b2))
#        print('c0:  %10.2e'%(c0))
        f1 = a0*tev**(-b0 -c0*np.log(tev))
        f2 = (1. + a2*tev**(-b2))/(1. +a1*tev**(-b1))
        RrRate = self.rrRate()
        Lev = em*ioneq*abund*f1*f2*tev*RrRate['rate']
        Lerg = em*ioneq*abund*f1*f2*const.boltzmann*temp*RrRate['rate']
        Lryd = em*ioneq*abund*f1*f2*(tev/const.ryd2Ev)*RrRate['rate']
        self.FreeBoundLossMao = {'f':f1*f2, 'LeV':Lev, 'Lerg':Lerg, 'Lryd':Lryd, 'rrRate':RrRate['rate'], 'temp':temp, 'tev':tev, 'abund':abund, 'ioneq':ioneq, 'em':em}

    def rrRate(self):
        """
        Provide the radiative recombination rate coefficient as a function of temperature (K).
        a revised copy of the Ion method
        """
        if hasattr(self, 'Temperature'):
            temperature = self.Temperature
        else:
            print(' temperature is not defined in rrRate')
            return

        rrparamsfile = util.ion2filename(self.IonStr) + '.rrparams'
        if hasattr(self, 'RrParams'):
            rrparams = self.RrParams
        elif os.path.isfile(rrparamsfile):
            self.RrParams = io.rrRead(self.IonStr)
            rrparams = self.RrParams
        else:
            self.RrRate = {'temperature':temperature, 'rate':np.zeros_like(temperature)}
            return

        if rrparams['rrtype'] == 1:
            a = rrparams['params'][3]
            b = rrparams['params'][4]
            t0 = rrparams['params'][5]
            t1 = rrparams['params'][6]
            rate = a/(np.sqrt(temperature/t0)*(1.+np.sqrt(temperature/t0))**(1.-b)*(1.+np.sqrt(temperature/t1))**(1.+b))
            return {'temperature':temperature, 'rate':rate}
        elif rrparams['rrtype'] == 2:
            a = rrparams['params'][3]
            b = rrparams['params'][4]
            t0 = rrparams['params'][5]
            t1 = rrparams['params'][6]
            c = rrparams['params'][7]
            t2 = rrparams['params'][8]
            b += c*np.exp(-t2/temperature)
            rate = a/(np.sqrt(temperature/t0)*(1.+np.sqrt(temperature/t0))**(1.-b)*(1.+np.sqrt(temperature/t1))**(1.+b))
            return {'temperature':temperature, 'rate':rate}
        elif rrparams['rrtype'] == 3:
            a = rrparams['params'][2]
            b = rrparams['params'][3]
            rate = a/(temperature/1.e+4)**b
            return {'temperature':temperature, 'rate':rate}
        else:
            return {'temperature':temperature, 'rate':np.zeros_like(temperature)}


    def mewe_gaunt_factor(self, **kwargs):
        """
        Calculate the Gaunt factor according to [106]_ for a single ion :math:`Z_z`.

        Using Eq. 9 of [106]_, the free-bound Gaunt factor for a single ion can be written as,

        .. math::
           G_{fb}^{Z,z} = \\frac{E_H}{k_BT}\\mathrm{Ab}(Z)\\frac{N(Z,z)}{N(Z)}f(Z,z,n)

        where :math:`E_H` is the ground-state potential of H, :math:`\mathrm{Ab}(Z)` is the
        elemental abundance, :math:`\\frac{N(Z,z)}{N(Z)}` is the fractional ionization, and
        :math:`f(Z,z,n)` is given by Eq. 10 and is approximated by Eq 16 as,

        .. math::
           f(Z,z,n) \\approx f_2(Z,z,n_0) = 0.9\\frac{\zeta_0z_0^4}{n_0^5}\exp{\left(\\frac{E_Hz_0^2}{n_0^2k_BT}\\right)} + 0.42\\frac{z^4}{n_0^{3/2}}\exp{\left(\\frac{E_Hz^2}{(n_0 + 1)^2k_BT}\\right)}

        where :math:`n_0` is the principal quantum number, :math:`z_0` is the effective charge (see Eq. 7 of [106]_),
        and :math:`\zeta_0` is the number of vacancies in the 0th shell and is given in Table 1 of [106]_.
        Here it is calculated in the same manner as in `fb_rad_loss.pro <http://www.chiantidatabase.org/idl/continuum/fb_rad_loss.pro>`_
        of the CHIANTI IDL library. Note that in the expression for :math:`G_{fb}`, we have not included
        the :math:`N_H/n_e` factor.

        Raises
        ------
        ValueError
            If no .fblvl file is available for this ion

        """
        # read in free-bound level information for the recombined ion
        # thermal energy scaled by H ionization potential
        scaled_energy = const.ryd2erg/const.boltzmann/self.Temperature
        # set variables used in Eq. 16 of Mewe et al.(1986)
        nameDict = util.convertName(self.IonStr)
        if not hasattr(self, 'Recombined_fblvl'):
            self.Recombined_fblvl = io.fblvlRead(nameDict['lower'])
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


    def freeBoundLoss(self,  includeAbund=True, includeIoneq=True, verner=True,  verbose=False):
        '''
        to calculate the free-bound (radiative recombination) energy loss rate coefficient of an ion,
        the ion is taken to be the target ion,
        including the elemental abundance and the ionization equilibrium population
        uses the Gaunt factors of [103]_ Karzas, W.J, Latter, R, 1961, ApJS, 6, 167
        provides rate = erg cm^-2 s^-1
        Notes
        -----
        - Uses the Gaunt factors of [103]_ for recombination to the ground level
        - Uses the photoionization cross sections of [2]_ to develop the free-bound cross section
        - Does not include the elemental abundance or ionization fraction
        - The specified ion is the target ion
        - uses the corrected version of the K-L bf gaunt factors available in CHIANTI V10
        - revised to calculate the bf cross, fb cross section and the maxwell energy distribution
        - the Verner cross sections are not included for now
        - using the RESTART formulation

        References
        ----------
        .. [2] Verner & Yakovlev, 1995, A&AS, `109, 125
            <http://adsabs.harvard.edu/abs/1995A%26AS..109..125V>`_
         '''
        if self.Defaults['wavelength'] != 'angstrom':
            print(' the continuum can only be calculated for wavelengths in angstroms')
            return

        xlabel =  self.Labels['radlossTlabel']
        ylabel = self.Labels['radlossYlabel']

        temperature = self.Temperature
#        tev = const.boltzmannEv*temperature
        em = self.Em

        if includeAbund:
            abund = self.Abundance
        else:
            abund = 1.
        if includeIoneq:
            ioneq = self.IoneqOne
        else:
            ioneq = np.ones_like(temperature)
        #
        # the target ion contains that data for fblvl
        #
        if hasattr(self,'Fblvl'):
            fblvl = self.Fblvl
            if 'errorMessage' in fblvl.keys():
                if verbose:
                    print(fblvl['errorMessage'])
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
                self.FreeBoundLoss = {'errorMessage':' no fblvl file for ion %s'%(self.IonStr)}
                if verbose:
                    print(' no fblvl file for ion %s'%(self.IonStr))
                return
        #
        #  need data for the recombined ion
        #
        if hasattr(self,'rFblvl'):
            rFblvl = self.rFblvl
        else:
            lower = self.nameDict['lower']
            lowerDict = util.convertName(lower)
            rFblvlname = lowerDict['filename'] +'.fblvl'
            if os.path.isfile(rFblvlname):
                self.rFblvl = io.fblvlRead(lower)
                rFblvl = self.rFblvl
            else:
                self.FreeBoundLoss = {'errorMessage':' no fblvl file for ion %s'%(lower)}
                if verbose:
                    print(' no fblvl file for ion %s'%(lower))
                return
        #
        # for the ionization potential, must use that of the recombined ion
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

        # get revised karzas-latter Gaunt factors
        klgbfn = chdata.Klgbfn
        #

        nTemp = temperature.size
        # constants for free bound
        # K_1 = 1.726 \times 10^{-28} = \frac{2^4 h e^2 }{3 \sqrt 3 m_e c}
        K1 = 2.**4*const.planck*const.q**2/(3.*np.sqrt(3.)*const.emass*const.light)
        # K_2 = 6.107 \times 10^5 = \frac{1}{ 2 m_e c^2 }
        K2 = 1./(2.*const.emass*const.light**2)

        #  \frac{1}{h c} (\frac{1}{2 m_e})^{1/2} (\frac{1}{\pi k})^{3/2}
        K3 = (1./(const.planck*const.light))*(1./(np.sqrt(2.*const.emass)))*(1./(np.pi*const.boltzmann)**1.5)
        K0 = 1.e-8*K1*K2*K3
#        K_4 = \frac{1}{\pi} (\frac{1}{2 \pi m_e k})^{1/2}
        K4 = (1./np.pi)*(1./(np.sqrt(2.*np.pi*const.emass*const.boltzmann)))

        fbLoss = np.zeros((nTemp),np.float64)
        ratg = np.zeros((nlvls),np.float64)
        iprLvlErg = np.zeros((nlvls),np.float64)
        gfGl = np.zeros((nTemp, nlvls, const.ngl), np.float64)
        gfInt = np.zeros((nTemp,  nlvls), np.float64)
        gfIntAllSum = np.zeros((nTemp), np.float64)
        egl = np.zeros((nTemp, nlvls, const.ngl), np.float64)
        scaledE = np.zeros((nTemp, nlvls, const.ngl), np.float64)

        peAll = []
        gfAll = []
        V1 = (1./(const.emass*const.light**2))*np.sqrt(2./(np.pi*const.emass*const.boltzmann))
        if verner:
            ilvl = 0
            lvl1 = 1
            ratg[0] = float(multr[0])/float(mult[0])
            iprLvlEv = self.Ipr - const.invCm2Ev*ecm[0]
            iprLvlErg[ilvl] = const.ev2Erg*iprLvlEv
            xgl = const.xgl
            wgl = const.wgl

            for itemp, atemp in enumerate(temperature):
                hnu = xgl*const.boltzmann*atemp + iprLvlErg[ilvl]
                wvl = 1.e+8*const.planck*const.light/hnu
                self.vernerCross(wvl)
                fbInt = wgl*(xgl*const.boltzmann*atemp + iprLvlErg[ilvl])**3*self.VernerCross
                fbLoss[itemp] = em[itemp]*abund*ioneq[itemp]*V1*ratg[0]*fbInt.sum()/np.sqrt(atemp)
        else:
            # without verner
            lvl1 = 0
        for ilvl in range(lvl1,nlvls):
            pqnIdx = pqn[ilvl] - 1
            lIdx = l[ilvl]
            klgbf = klgbfn[pqnIdx]
            pe = klgbf['pe']
            peAll.append(pe)
            gf = klgbf['klgbf'][lIdx]
            gfAll.append(gf)
            iprLvlEv = self.Ipr - const.invCm2Ev*ecm[ilvl]
            iprLvlErg[ilvl] = const.ev2Erg*iprLvlEv
            ratg[ilvl] = float(multr[ilvl])/float(mult[0]) # ratio of statistical weights

            for itemp,  atemp in enumerate(temperature):

                egl[itemp, ilvl] = const.xgl*const.boltzmannEv*atemp +  iprLvlEv
                scaledE[itemp, ilvl] = egl[itemp,  ilvl]/self.Ipr
                tck = splrep(np.log(pe), np.log(gf),  s=0)
                gflog = splev(np.log(scaledE[itemp, ilvl]), tck,  der=0,  ext=1)
                gfGl[itemp, ilvl] = np.exp(gflog)
                for igl, awgl in enumerate(const.wgl):
                    gfInt[itemp, ilvl] += awgl*gfGl[itemp, ilvl, igl]
                gfIntSum = gfInt.sum(axis=1)
        for itemp,  atemp in enumerate(temperature):
            for ilvl in range(lvl1,nlvls):
                gfIntAllSum[itemp] += ratg[ilvl]*iprLvlErg[ilvl]**2*gfIntSum[itemp]/float(pqn[ilvl])
        fbkl = em*ioneq*abund*K1*K2*K4*gfIntAllSum/np.sqrt(temperature)
        fbLoss += fbkl

        self.FreeBoundLoss = {'rate':fbLoss, 'fbkl':fbkl, 'temperature':temperature,
            'xlabel':xlabel, 'ylabel':ylabel}


    def freeBoundRate(self, verner=True, verbose=False):
        """
        Calculates the free-bound (radiative recombination) rate of an ion.
        Provides emissivity in units of ergs :math:`\mathrm{cm}^{-2}` :math:`\mathrm{s}^{-1}` :math:`\mathrm{str}^{-1}` :math:`\mathrm{\AA}^{-1}` for an individual ion.  If includeAbund is set,
        the abundance is included.  If includeIoneq is set, the ionization equililibrium for the given
        ion is included
        This method is only for the purposes of testing.
        The best rr/FB rates come from the ion class

        Notes
        -----
        - Uses the Gaunt factors of [103]_ for recombination to the excited levels
        - Uses the photoionization cross sections of [3]_ to develop the free-bound cross section
        - Include the elemental abundance and ionization fraction by default
        - The specified ion is the target ion
        - uses the corrected version of the K-L bf gaunt factors available in CHIANTI V10
        - revised to calculate the bf cross, fb cross section and the maxwell energy distribution

        """

        if self.Defaults['wavelength'] != 'angstrom':
            print(' the continuum can only be calculated for wavelengths in angstroms')
            return

        temperature = self.Temperature

        ngl = 15
        zgl = special.roots_laguerre(ngl)
        xgl = zgl[0]
        wgl = zgl[1]

        xlabel = self.Labels['rateTlabel']
        ylabel = self.Labels['rateYlabel']

#        if np.array_equal(self.Em, np.ones_like(self.Em)):
#            ylabel += '($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'


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
            self.Fblvl = {'mult':[1., 1.]}
            fblvl = self.Fblvl
        else:
            fblvlname = self.nameDict['filename']+'.fblvl'
            if os.path.isfile(fblvlname):
                self.Fblvl = io.fblvlRead(self.IonStr)
                fblvl = self.Fblvl
            # in case there is no fblvl file
            else:
                if verbose:
                    print(' No fblvl file for ion %s'%(self.IonStr))
                self.FreeBound = {'errorMessage':' No fblvl file for ion %s'%(self.IonStr)}
                return
        #
        #  need data for the recombined ion
        #
        if hasattr(self,'rFblvl'):
            rFblvl = self.rFblvl
        else:
            lower = self.nameDict['lower']
            lowerDict = util.convertName(lower)
            rFblvlname = lowerDict['filename'] +'.fblvl'
            if os.path.isfile(rFblvlname):
                self.rFblvl = io.fblvlRead(lower)
                rFblvl = self.rFblvl
            else:
                if verbose:
                    print(' No fblvl file for ion %s'%(lower))
                self.FreeBound = {'errorMessage':' no fblvl file for ion %s'%(lower)}
                return
        #
        # for the ionization potential, must use that of the recombined ion
        #
        nlvls = len(rFblvl['lvl'])
        # pqn = principle quantum no. n
        pqn = np.asarray(rFblvl['pqn'], 'int64')
        conf = rFblvl['conf']

        # l is angular moment quantum no. L
        l = rFblvl['l']
        # energy level in inverse cm
        ecm = rFblvl['ecm']
        # statistical weigths/multiplicities
        multr = rFblvl['mult']
        mult = fblvl['mult']

        # get revised karzas-latter Gaunt factors
        klgbfn = chdata.Klgbfn
        #
        nWvl = ngl
        nTemp = temperature.size
        #
        #  the verner cross-section is included in this version
#        if verner:
#            self.vernerCross(wvl)
#            vCross = self.VernerCross
#            lvl1 = 1
#        else:
#            lvl1 = 0
#        lvl1 = 0
#        expfun = np.zeros((nlvls,nTemp,nWvl),np.float64)
#        mask = np.zeros((nlvls,nTemp,nWvl),np.bool_)
        fbn = np.zeros((nlvls,nTemp,nWvl),np.float64)
        fbIntensity = np.zeros((nlvls,nTemp,nWvl),np.float64)
        ratg = np.zeros((nlvls),np.float64)
#        mygf = np.zeros((nlvls), np.float64)
#        ratg[0] = float(multr[0])/float(mult[0])
        iprLvlEv = self.Ipr - const.invCm2Ev*ecm[0]
        iprLvlErg = const.ev2Erg*iprLvlEv
        edgeLvlAng = []
        mygfInt = np.zeros(nlvls,  np.float64)
        fbRate = np.zeros((nlvls,  nTemp), np.float64)
#        hnu = const.planck*const.light/(1.e-8*wvl)
#        hnuEv = const.ev2Ang/wvl
        for ilvl in range(nlvls):
            ratg[ilvl] = float(multr[ilvl])/float(mult[0]) # ratio of statistical weights

        # constants for free bound
        # K_1 = 1.726 \times 10^{-28} = \frac{2^4 h e^2 }{3 \sqrt 3 m_e c}
        K1 = 2.**4*const.planck*const.q**2/(3.*np.sqrt(3.)*const.emass*const.light)
        # K_2 = 6.107 \times 10^5 = \frac{1}{ 2 m_e c^2 }
        K2 = 1./(2.*const.emass*const.light**2)

        #  \frac{1}{h c} (\frac{1}{2 m_e})^{1/2} (\frac{1}{\pi k})^{3/2}
#        K3 = (1./(const.planck*const.light))*(1./(np.sqrt(2.*const.emass)))*(
#            1./(np.pi*const.boltzmann)**1.5)
        K5 = np.sqrt(8./(np.pi*const.emass*const.boltzmann))
        K0 = K1*K2*K5

        # restrict calculations to temperatures where the ionization equ. is non-zero
#        goodT = [ it for it, anioneq in enumerate(self.IoneqOne) if anioneq > 0.]
        if verner:
#            vrrRate = np.zeros((nTemp), np.float64)
            V3 = (4.*np.pi)/(const.light**2)*(1./(2.*np.pi*const.emass))**(1.5)*np.sqrt(1./const.boltzmann)
            ilvl = 0
            lvl1 = 1
            iprLvlEv = self.Ipr - const.invCm2Ev*ecm[ilvl]
            edgeLvlAng.append(const.ev2Ang/iprLvlEv)

            for it, atemp in enumerate(temperature):
                xkt = xgl*const.boltzmann*atemp
                hnu = xkt + iprLvlErg
                wvl = (1.e8*const.planck*const.light)/hnu
                self.vernerCross(wvl)
                vernerSum = (wgl*(hnu)**2*self.VernerCross).sum()
                fbRate[0, it] = V3*ratg[0]*vernerSum/np.sqrt(atemp)
#                print(' verner %5i %5i  %10.2e'%(ilvl, itemp, atemp))
#                mask[0,itemp] = 1.e+8/wvl < (self.IprCm - ecm[0])
#                expfun[0, itemp] = np.exp((iprLvlErg - 1.e+8*const.planck*const.light/wvl)/
#                    (const.boltzmann*atemp))
#                fbn[0, itemp] = (const.planck*const.light/(1.e-8*wvl))**5 \
#                    *const.verner*ratg[0]*expfun[0,itemp]*vCross/atemp**1.5
#                fbIntensity[0,itemp] = abund*ioneq[itemp]*em[itemp]*fbn[0,  itemp]

        else:
            lvl1 = 0

        for ilvl in range(lvl1,nlvls):
            pqnIdx = pqn[ilvl] - 1
            lIdx = l[ilvl]
            klgbf = klgbfn[pqnIdx]
            pe = klgbf['pe']
            gf = klgbf['klgbf'][lIdx]
            iprLvlEv = self.Ipr - const.invCm2Ev*ecm[ilvl]
            edgeLvlAng.append(const.ev2Ang/iprLvlEv)
            iprLvlErg = const.ev2Erg*iprLvlEv

            # scaled energy is relative to the ionization potential for each level

#            scaledE = hnuEv/iprLvlEv
#            badLong = hnuEv < iprLvlEv
#            if verbose:
#                print(' %i nWvl  %i badlong.sum %i'%(ilvl, nWvl,  badLong.sum()))
#            if badLong.sum() < nWvl:
            tck = splrep(np.log(pe), np.log(gf),  s=0)
#            gflog = splev(np.log(scaledE), tck,  der=0,  ext=3)
#            mygf[ilvl] = np.exp(gflog)
#            mygf[ilvl][badLong] = 0.
#                if verbose:
#                    print(' ilvl %i   mygf.max() %10.2e'%(ilvl, np.max(mygf[ilvl])))
#            ratg[ilvl] = float(multr[ilvl])/float(mult[0]) # ratio of statistical weights

            for itemp,  atemp in enumerate(self.Temperature):
                xkt = xgl*const.boltzmann*atemp
                hnu = xkt + iprLvlErg
                hnuEv = hnu/const.ev2Erg
                scaledE = hnuEv/iprLvlEv
#                badLong = hnuEv < iprLvlEv
#                tck = splrep(np.log(pe), np.log(gf),  s=0)
                gflog = splev(np.log(scaledE), tck,  der=0,  ext=3)
                mygf = np.exp(gflog)
                mygfInt = (wgl*mygf/hnu).sum()
                fbRate[ilvl,  itemp] = (K0*ratg[ilvl]*iprLvlErg**2*mygfInt/float(pqn[ilvl]))/np.sqrt(atemp)
#                    iprLvlErg - hnu)/(const.boltzmann*atemp
#                    print(' kl %5i %5i  %10.2e'%(ilvl, itemp, atemp))
                # this test if not really needed
#                    if ioneq[itemp] > 0.:
#                    expfun = np.exp((iprLvlErg - hnu)/(const.boltzmann*atemp))
#                    expfunIsnan = np.isnan(expfun)
#                    expfunIsinf = np.isinf(expfun)
#                    if expfunIsnan.sum() > 0:
#                        print('%i some of expfun are nan %i'%(itemp, expfunIsnan.sum()))
#                        expfun[expfunIsnan] = 0.
#                    elif expfunIsinf.sum() > 0:
#                        print('%i some of expfun are inf %i'%(itemp, expfunIsinf.sum()))
#                        expfun[expfunIsinf] = 0.

#                if self.Defaults['flux'] == 'energy':
#                    fbn[ilvl, itemp] = K0*hnu**2*expfun*iprLvlErg**2*ratg[ilvl]*mygf[ilvl]  \
#                        /(atemp**(1.5)*float(pqn[ilvl]))
#                elif self.Defaults['flux'] == 'photon':
#                    fbn[ilvl, itemp] = K0*hnu*expfun*iprLvlErg**2*ratg[ilvl]*mygf[ilvl]  \
#                        /(atemp**(1.5)*float(pqn[ilvl]))

#                fbnIsNan = np.isnan(fbn[ilvl, itemp])
#                fbnIsInf = np.isinf(fbn[ilvl, itemp])

#                if fbnIsNan.sum() > 0:
#                    print('%s  %i some fbn are Nan %i'%(self.IonStr, itemp, fbnIsNan.sum()))
#                    fbn[ilvl, itemp,  fbnIsNan] = 0.
#                elif fbnIsInf.sum() > 0:
#                    print('%s %i some fbn are Inf %i'%(self.IonStr, itemp, fbnIsInf.sum()))
#                    fbn[ilvl, itemp,  fbnIsInf] = 0.
#                mask[ilvl,itemp] = 1.e+8/wvl < (self.IprCm - ecm[ilvl])
#                fbIntensity[ilvl, itemp] = abund*em[itemp]*ioneq[itemp]*fbn[ilvl, itemp]
#
        fbTest = fbRate.sum(axis=0)
        fbTestisNan = np.isnan(fbTest)
        fbTestisInf = np.isinf(fbTest)
        if fbTestisNan.sum() > 0:
            print('fb has some Nans %i'%(fbTestisNan.sum()))
        if fbTestisInf.sum() > 0:
            print('fb has some Infs %i'%(fbTestisInf.sum()))
#        fb = np.ma.array(fbIntensity.sum(axis=0))
        fb = np.ma.array(fbRate.sum(axis=0))
#        fb_mask = fb <= 0. or fb == np.inf
        fb_mask = np.where(fb.data is np.nan,  1,  0)
        fb.mask = fb_mask
        fb.fill_value = 0.
        #
#        self.FreeBoundRate = {'rate':fbRate, 'temperature':temperature, 'gf':mygf, 'edgeLvlAng':edgeLvlAng,  'fbn':fbn,
#            'xlabel':xlabel, 'ylabel':ylabel,  'pqn':pqn, 'conf':conf}
#
        self.FreeBoundRate = {'rate':fb, 'fbLvl': fbRate,  'temperature':temperature, 'gf':mygf, 'edgeLvlAng':edgeLvlAng,  'fbn':fbn,
            'xlabel':xlabel, 'ylabel':ylabel,  'pqn':pqn, 'conf':conf}



    def freeBound(self, wvl, includeAbund = True, includeIoneq = True, verner=True, verbose=False):
        """
        Calculates the free-bound (radiative recombination) continuum emissivity of an ion.
        Provides emissivity in units of ergs :math:`\mathrm{cm}^{-2}` :math:`\mathrm{s}^{-1}` :math:`\mathrm{str}^{-1}` :math:`\mathrm{\AA}^{-1}` for an individual ion.  If includeAbund is set,
        the abundance is included.  If includeIoneq is set, the ionization equililibrium for the given
        ion is included

        Notes
        -----
        - Uses the Gaunt factors of [103]_ for recombination to the excited levels
        - Uses the photoionization cross sections of [3]_ to develop the free-bound cross section
        - Include the elemental abundance and ionization fraction by default
        - The specified ion is the target ion
        - uses the corrected version of the K-L bf gaunt factors available in CHIANTI V10
        - revised to calculate the bf cross, fb cross section and the maxwell energy distribution

        References
        ----------
        .. [3] Verner & Yakovlev, 1995, A&AS, `109, 125
            <http://adsabs.harvard.edu/abs/1995A%26AS..109..125V>`_
        """

        if self.Defaults['wavelength'] != 'angstrom':
            print(' the continuum can only be calculated for wavelengths in angstroms')
            return

        temperature = self.Temperature
        wvl = np.asarray(wvl, np.float64)
        em = self.Em

        nWvl = wvl.size
        nTemp = temperature.size


        xlabel = self.Labels['xlabel']
        ylabel = self.Labels['spectrumYlabel']

        if np.array_equal(self.Em, np.ones_like(self.Em)):
            ylabel += '($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'

        if includeAbund:
            abund = self.Abundance
        else:
            abund = 1.

        if includeIoneq:
            ioneq = self.IoneqOne
            goodT = [it for it in range(self.Ntemp) if self.IoneqOne[it] > 0.]
        else:
            ioneq = np.ones_like(temperature)
            goodT = range(nTemp)
        #
        # the target ion contains the data for fblvl
        #
        if hasattr(self,'Fblvl'):
            fblvl = self.Fblvl
            if 'errorMessage' in fblvl.keys():
                self.FreeBound = fblvl
                return
        elif self.Z == self.Stage-1:
            #dealing with the fully ionized stage
            self.Fblvl = {'mult':[1., 1.]}
            fblvl = self.Fblvl
        else:
            fblvlname = self.nameDict['filename']+'.fblvl'
            if os.path.isfile(fblvlname):
                self.Fblvl = io.fblvlRead(self.IonStr)
                fblvl = self.Fblvl
            # in case there is no fblvl file
            else:
                if verbose:
                    print(' No fblvl file for ion %s'%(self.IonStr))
                self.FreeBound = {'errorMessage':' No fblvl file for ion %s'%(self.IonStr)}
                return
        #
        #  need data for the recombined ion
        #
        if hasattr(self,'rFblvl'):
            rFblvl = self.rFblvl
        else:
            lower = self.nameDict['lower']
            lowerDict = util.convertName(lower)
            rFblvlname = lowerDict['filename'] +'.fblvl'
            if os.path.isfile(rFblvlname):
                self.rFblvl = io.fblvlRead(lower)
                rFblvl = self.rFblvl
            else:
                if verbose:
                    print(' No fblvl file for ion %s'%(lower))
                self.FreeBound = {'errorMessage':' no fblvl file for ion %s'%(lower)}
                return
        #
        # for the ionization potential, must use that of the recombined ion
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

        # get revised karzas-latter Gaunt factors
        klgbfn = chdata.Klgbfn
        #
        #

        expfun = np.zeros((nlvls,nTemp,nWvl),np.float64)
        fbn = np.zeros((nlvls,nTemp,nWvl),np.float64)
        fbIntensity = np.zeros((nlvls,nTemp,nWvl),np.float64)
        ratg = np.zeros((nlvls),np.float64)
        mygf = np.zeros((nlvls, nWvl))
        ratg[0] = float(multr[0])/float(mult[0])
        iprLvlEv = self.Ipr - const.invCm2Ev*ecm[0]
        iprLvlErg = const.ev2Erg*iprLvlEv
        edgeLvlAng = []

        hnu = const.planck*const.light/(1.e-8*wvl)
        hnuEv = const.ev2Ang/wvl

        # constants for free bound
        # K_1 = 1.726 \times 10^{-28} = \frac{2^4 h e^2 }{3 \sqrt 3 m_e c}
        K1 = 2.**4*const.planck*const.q**2/(3.*np.sqrt(3.)*const.emass*const.light)
        # K_2 = 6.107 \times 10^5 = \frac{1}{ 2 m_e c^2 }
        K2 = 1./(2.*const.emass*const.light**2)

        #  \frac{1}{h c} (\frac{1}{2 m_e})^{1/2} (\frac{1}{\pi k})^{3/2}
        K3 = (1./(const.planck*const.light))*(1./(np.sqrt(2.*const.emass)))*(1./(np.pi*const.boltzmann)**1.5)
        K0 = 1.e-8*K1*K2*K3

#        V1 = np.sqrt(2.)*const.planck**4*const.light**2/(4*(np.pi*const.emass*const.boltzmann)**(1.5))
        # restrict calculations to temperatures where the ionization equ. is non-zero


        #  the verner yakovlev photionizaton cross sections

        if verner:
            lvl1 = 1
            self.vernerCross(wvl)
            ilvl = 0
            iprLvlEv = self.Ipr - const.invCm2Ev*ecm[ilvl]
            edgeLvlAng.append(const.ev2Ang/iprLvlEv)

            for itemp in  goodT:
                atemp = temperature[itemp]

                xponent = (iprLvlErg - hnu)/(const.boltzmann*atemp)
                expfun[0, itemp] = np.where(xponent <= 0.,  np.exp(xponent),  0.)

                c1 = const.verner*ratg[0]*expfun[0,itemp]*self.VernerCross/atemp**1.5

                if self.Defaults['flux'] == 'energy':
                    fbn[0,itemp] = (const.planck*const.light/(1.e-8*wvl))**5*c1
                elif self.Defaults['flux'] == 'photon':
                    fbn[0,itemp] = (const.planck*const.light/(1.e-8*wvl))**4*c1

                fbIntensity[ilvl, itemp] = em[itemp]*fbn[ilvl, itemp]

        else:
            lvl1 = 0

        # the karzas latter gaunt factors are used

        for ilvl in range(lvl1,nlvls):
            pqnIdx = pqn[ilvl] - 1
            lIdx = l[ilvl]
            klgbf = klgbfn[pqnIdx]
            pe = klgbf['pe']
            gf = klgbf['klgbf'][lIdx]
            iprLvlEv = self.Ipr - const.invCm2Ev*ecm[ilvl]
            edgeLvlAng.append(const.ev2Ang/iprLvlEv)
            iprLvlErg = const.ev2Erg*iprLvlEv

            # scaled energy is relative to the ionization potential

            scaledE = hnuEv/self.Ipr

            tck = splrep(np.log(pe), np.log(gf),  s=0)
            gflog = splev(np.log(scaledE), tck,  der=0,  ext=3)
            mygf[ilvl] = np.where(hnuEv >= iprLvlEv, np.exp(gflog),  0.)

            ratg[ilvl] = float(multr[ilvl])/float(mult[0]) # ratio of statistical weights

            for itemp in  goodT:
                atemp = temperature[itemp]

                xponent = (iprLvlErg - hnu)/(const.boltzmann*atemp)
                expfun = np.where(xponent <= 0.,  np.exp(xponent),  0.)

                if self.Defaults['flux'] == 'energy':
                    fbn[ilvl, itemp] = K0*hnu**2*expfun*iprLvlErg**2*ratg[ilvl]*mygf[ilvl]  \
                        /(atemp**(1.5)*float(pqn[ilvl]))
                elif self.Defaults['flux'] == 'photon':
                    fbn[ilvl, itemp] = K0*hnu*expfun*iprLvlErg**2*ratg[ilvl]*mygf[ilvl]  \
                        /(atemp**(1.5)*float(pqn[ilvl]))


                fbIntensity[ilvl, itemp] = em[itemp]*fbn[ilvl, itemp]

        fb = fbIntensity.sum(axis=0)

        if includeAbund:
            fb *= abund

        if includeIoneq:
            for itemp, one in enumerate(self.IoneqOne):
                fb[itemp] *= one
        #
        self.FreeBound = {'intensity':fb.squeeze(), 'temperature':temperature,'wvl':wvl, 'em':em, \
            'abund':abund, 'ioneq':ioneq, 'gf':mygf, 'edgeLvlAng':edgeLvlAng,  'fbn':fbn.squeeze(),
            'xlabel':xlabel, 'ylabel':ylabel}


    def vernerCross(self, wvl):
        """
        Calculates the photoionization cross-section using data from [102]_ for
        transitions to the ground state.

        The photoionization cross-section can be expressed as :math:`\sigma_i^{fb}=F(E/E_0)` where
        :math:`F` is an analytic fitting formula given by Eq. 1 of [102]_,

        .. math::
           F(y) = ((y-1)^2 + y_w^2)y^{-Q}(1 + \sqrt{y/y_a})^{-P},

        where :math:`E` is the photon energy, :math:`n` is the principal quantum number,
        :math:`l` is the orbital quantum number, :math:`Q = 5.5 + l - 0.5P`, and
        :math:`\sigma_0,E_0,y_w,y_a,P` are fitting paramters. These can be read in using
        `ChiantiPy.tools.io.vernerRead`.

        Parameters:

        wvl : `list`, `ndarray`
            wavelength in Angstroms

        """
        # read verner data
        verner_info = chdata.Verner_info
        eth = verner_info['eth'][self.Z,self.Stage-1]   #*const.ev2Erg
        yw = verner_info['yw'][self.Z,self.Stage-1]
        ya = verner_info['ya'][self.Z,self.Stage-1]
        p = verner_info['p'][self.Z,self.Stage-1]

        # convert from megabarn to cm^2
        sigma0 = verner_info['sig0'][self.Z,self.Stage-1]*1e-18
        e0 = verner_info['e0'][self.Z,self.Stage-1]  #*const.ev2Erg
        q = 5.5 + verner_info['l'][self.Z,self.Stage-1] - 0.5*p

        # scaled photon energy
        en = const.ev2Ang/np.asarray(wvl)
        y = en/e0
        # fitting function
        F = ((y - 1.)**2 + yw**2)*(y**(-q))*(1. + np.sqrt(y/ya))**(-p)
        cross_section = sigma0*F

        self.VernerCross = np.where(en < eth, 0., cross_section)

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
