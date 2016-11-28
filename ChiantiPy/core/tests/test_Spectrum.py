"""
Tests for the spectrum and bunch classes
"""

import numpy as np
import pytest

from ChiantiPy.core import spectrum,bunch

#set temperature, density, wavelength
temperature_1 = np.array([1e+6,4e+6,1e+7])
temperature_2 = np.logspace(5,8,10)
density = 1e+9
wavelength = np.linspace(10,100,1000)
wavelength_range = [wavelength[0],wavelength[-1]]
min_abund=1.e-4
ion_list = ['fe_15','fe_16']

# test the spectrum class
def test_spectrum():
    _tmp_spec = spectrum(temperature_1,density,wavelength,minAbund=min_abund)
    # TODO: need to assert something here, not clear what exactly to test yet

# test the bunch class
def test_bunch():
    _tmp_bunch = bunch(temperature_2,density,wvlRange=wavelength_range,ionList=ion_list)
    # TODO: need to assert something here, not clear what exactly to test yet
