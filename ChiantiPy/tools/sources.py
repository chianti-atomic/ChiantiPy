"""
Blackbody temperature calculations
"""

import numpy as np

import ChiantiPy.tools.constants as const

class blackStar:
    """
    Calculate blackbody radiation

    Parameters
    ----------
    temperature : `~numpy.ndarray`
        Temperature in Kelvin
    radius : `~numpy.ndarray`
        Stellar radius in cm

    Attributes
    ----------
    Temperature : `~numpy.ndarray`
        Temperature in Kelvin
    Radius : `~numpy.ndarray`
        Stellar radius in cm
    Incident : `~numpy.ndarray`
        Blackbody photon distribution
    """

    def __init__(self, temperature, radius):
        self.Temperature = temperature
        self.Radius = radius

    def incident(self, distance, energy):
        """
        Calculate photon distribution times the visible cross-sectional area.

        Parameters
        ----------
        distance : `~numpy.ndarray`
            Distance to the stellar object in cm
        energy : `~numpy.ndarray`
            Energy range in erg

        Notes
        -----
        This function returns the photon distribution instead of the distribution times the cross-sectional area. Is this correct? Why is the incident photon distribution calculated at all?
        """
        print((' distance %10.2e  energy '%(energy)))
        bb = blackbody(self.Temperature, energy)
        out = const.pi*(self.Radius/distance)**2*bb['photons']
        self.Incident = bb


def blackbody(temperature, variable, hnu=1):
    """
    Calculate the blackbody photon distribution as a function of energy (`hnu` = 1) or as a function of wavelength (`hnu` = 0) in units of  :math:`\mathrm{photons}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}\,\mathrm{str}^{-1}\,\mathrm{erg}^{-1}`

    Parameters
    ----------
    temperature : `~numpy.float64`
        Temperature at which to calculate the blackbody photon distribution
    variable : `~numpy.ndarray`
        Either energy (in erg) or wavelength (in angstrom)
    hnu : `int`
        If 1, calculate distribution as a function of energy. Otherwise, calculate it as a function of wavelength

    Returns
    -------
    {'photons', 'temperature', 'energy'} or {'photons', 'temperature', 'wvl'} : `dict`
    """
    if hnu:
        energy = variable
        bb =(2./(const.planck*(const.hc**2)))*energy**2/(np.exp(energy/(const.boltzmann*temperature)) - 1.)
        return {'photons':bb, 'temperature':temperature, 'energy':energy}
    else:
        wvl = 1.e-8*variable
        bb = ((2.*const.pi*const.light)/wvl**4)/(np.exp(const.hc/(wvl*const.boltzmann*temperature)) - 1.)
        return {'photons':bb, 'temperature':temperature, 'wvl':wvl}
