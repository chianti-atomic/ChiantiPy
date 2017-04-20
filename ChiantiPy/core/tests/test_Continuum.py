"""
Tests for the continuum class
"""

import numpy as np
import pytest

from ChiantiPy.core import Continuum

# test ion
test_ion = 'fe_15'
# set some wavelength and temperature arrays
temperature_1 = 1e+6
temperature_2 = np.logspace(5,8,10)
wavelength = np.linspace(10,100,100)
# create continuum object for testing
tmp_cont = Continuum(test_ion,temperature_2)
# create continuum object for which there is no free-bound information available
tmp_cont_no_fb = Continuum('fe_3',temperature_2)

# test temperature input
def test_temperature():
    _tmp_cont = Continuum(test_ion,temperature_1)
    assert np.array(temperature_1) == _tmp_cont.temperature
    _tmp_cont = Continuum(test_ion,temperature_2)
    assert np.all(temperature_2 == _tmp_cont.temperature)

# test the free-free calculation
def test_free_free():
    # call free-free emission and loss rate methods
    tmp_cont.calculate_free_free_emission(wavelength)
    assert hasattr(tmp_cont,'free_free_emission')
    tmp_cont.calculate_free_free_loss()
    assert hasattr(tmp_cont,'free_free_loss')

# test the free-bound calculation
def test_free_bound():
    # free-bound emission and loss rate methods
    tmp_cont.calculate_free_bound_emission(wavelength)
    assert hasattr(tmp_cont,'free_bound_emission')
    tmp_cont.calculate_free_bound_loss()
    assert hasattr(tmp_cont,'free_bound_loss')
    # test an error being raised if no free-bound information is available
    with pytest.raises(ValueError, message='Expecting ValueError when no free-bound information is available'):
        tmp_cont_no_fb.calculate_free_bound_emission(wavelength)
