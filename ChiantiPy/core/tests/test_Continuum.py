"""
Tests for the continuum class
"""

import numpy as np
import pytest

from ChiantiPy.core import continuum

# test ion
test_ion = 'fe_15'
# set some wavelength and temperature arrays
temperature_1 = 1e+6
temperature_2 = np.logspace(5,8,10)
wavelength = np.linspace(10,100,100)
# create continuum object for testing
tmp_cont = continuum(test_ion,temperature_2)


# test temperature input
def test_temperature():
    _tmp_cont = continuum(test_ion,temperature_1)
    assert np.array(temperature_1)==temperature_1
    _tmp_cont = continuum(test_ion,temperature_2)
    assert np.all(temperature_2==tmp_cont.Temperature)

# test the free-free calculation
def test_freeFree():
    # call free-free continuum method
    tmp_cont.freeFree(wavelength)
    assert hasattr(tmp_cont,'FreeFree')
    # check that all expected keys are there
    # FIXME: need to check 'ff' key too but conditions when it should be there are not clear
    assert all((k in tmp_cont.FreeFree for k in ('rate','temperature','wvl') ))

# test the free-bound calculation
def test_freeBound():
    # free-bound continuum methods
    tmp_cont.freeBound(wavelength)
    assert hasattr(tmp_cont,'FreeBound')
    # check that all expected keys are there
    assert all((k in tmp_cont.FreeBound for k in ('rate','temperature','wvl') ))
