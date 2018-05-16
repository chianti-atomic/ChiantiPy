"""
Tests for the continuum class
"""

import numpy as np
import pytest

from ChiantiPy.core import Continuum

# test ion
test_ion = 'fe_15'
# set some wavelength and temperature arrays
temperature_scalar = 1e+6
temperature_array = np.logspace(5,8,10)
wavelength_array = np.linspace(10,100,100)
# create continuum object for testing
tmp_cont_scalar = Continuum.continuum(test_ion, temperature_scalar)
tmp_cont_array = Continuum.continuum(test_ion, temperature_array)
# create continuum object for which there is no free-bound information available
tmp_cont_no_fb = Continuum.continuum('fe_3', temperature_array)


def test_temperature():
    _tmp_cont = Continuum.continuum(test_ion, temperature_scalar)
    assert np.all(np.atleast_1d(temperature_scalar) == _tmp_cont.Temperature)
    _tmp_cont = Continuum.continuum(test_ion, temperature_array)
    assert np.all(temperature_array == _tmp_cont.Temperature)


def test_free_free_scalar():
    # call free-free emission and loss rate methods for scalar temperature
    # emission--wavelength array
    tmp_cont_scalar.freeFree(wavelength_array)
    assert hasattr(tmp_cont_scalar, 'FreeFree')
    assert tmp_cont_scalar.FreeFree['intensity'].shape == wavelength_array.shape
    # loss
    tmp_cont_scalar.freeFreeLoss()
    assert hasattr(tmp_cont_scalar, 'FreeFreeLoss')
    assert tmp_cont_scalar.FreeFreeLoss['rate'].shape == (1,)


def test_free_free_array():
    # call free-free emission and loss rate methods for scalar temperature
    # emission--wavelength array
    tmp_cont_array.freeFree(wavelength_array)
    assert hasattr(tmp_cont_array, 'FreeFree')
    assert tmp_cont_array.FreeFree['intensity'].shape == temperature_array.shape+wavelength_array.shape
    # loss
    tmp_cont_array.freeFreeLoss()
    assert hasattr(tmp_cont_array, 'FreeFreeLoss')
    assert tmp_cont_array.FreeFreeLoss['rate'].shape == temperature_array.shape


def test_free_bound_scalar():
    # free-bound emission and loss rate methods for scalar temperature
    # emiss--wavelength array
    tmp_cont_scalar.freeBound(wavelength_array)
    assert hasattr(tmp_cont_scalar,'FreeBound')
    assert tmp_cont_scalar.FreeBound['intensity'].shape == wavelength_array.shape
    # loss
    tmp_cont_scalar.freeBoundLoss()
    assert hasattr(tmp_cont_scalar, 'FreeBoundLoss')
    assert tmp_cont_scalar.FreeBoundLoss['rate'].shape == (1,)


def test_free_bound_array():
    # free-bound emission and loss rate methods for temperature array
     # emiss--wavelength array
    tmp_cont_array.freeBound(wavelength_array)
    assert hasattr(tmp_cont_array,'FreeBound')
    assert tmp_cont_array.FreeBound['intensity'].shape == temperature_array.shape+wavelength_array.shape
    # loss
    tmp_cont_array.freeBoundLoss()
    assert hasattr(tmp_cont_array, 'FreeBoundLoss')
    assert tmp_cont_array.FreeBoundLoss['rate'].shape == temperature_array.shape


def test_free_bound_no_info():
    # raise error if no free-bound information is available
    tmp_cont_no_fb.freeBound(wavelength_array)
    assert 'errorMessage' in tmp_cont_no_fb.FreeBound.keys()
