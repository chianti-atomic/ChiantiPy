import numpy as np
import chianti.constants as const
class blackStar:
    '''temperature in K, radius is the stellar radius in cm '''
    def __init__(self, temperature, radius):
        self.Temperature = temperature
        self.Radius = radius
    def incident(self, distance, energy):
        ''' distance in cm and energy in erg'''
        print((' distance %10.2e  energy '%(energy)))
        bb = blackbody(self.Temperature, energy)
        out = const.pi*(self.Radius/distance)**2*bb['photons']
        self.Incident = bb
    #
    # ---------------------------------------------------------------------
    #
def blackbody(temperature, variable, hnu=1):
    ''' to calculate the black body photon distribution as a function of energy in erg (hnu = 1) or as a function
    of wavelength (Angstroms) (hnu=0)
    photons cm^-2 s^-1 str^-1 ergs^-1'''
    if hnu:
        energy = variable
        bb =(2./(const.planck*(const.hc**2)))*energy**2/(np.exp(energy/(const.boltzmann*temperature)) - 1.)
        return {'photons':bb, 'temperature':temperature, 'energy':energy}
    else:
        wvl = 1.e-8*variable
        bb = ((2.*const.pi*const.light)/wvl**4)/(np.exp(const.hc/(wvl*const.boltzmann*temperature)) - 1.)
        return {'photons':bb, 'temperature':temperature, 'wvl':wvl}
