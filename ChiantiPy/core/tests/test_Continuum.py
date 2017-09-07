"""
Tests for the continuum class
"""

import numpy as np
import pytest

from ChiantiPy.core import continuum

# test ion
test_ion = 'fe_15'
# set some wavelength and temperature arrays
temperature_scalar = 1e+6
temperature_array = np.logspace(5,8,10)
wavelength_scalar = 50.0
wavelength_array = np.linspace(10,100,100)
# create continuum object for testing
tmp_cont_scalar = continuum(test_ion, temperature_scalar)
tmp_cont_array = continuum(test_ion, temperature_array)
# create continuum object for which there is no free-bound information available
tmp_cont_no_fb = continuum('fe_3', temperature_array)


def test_temperature():
    _tmp_cont = Continuum(test_ion, temperature_scalar)
    assert np.all(np.atleast_1d(temperature_scalar) == _tmp_cont.Temperature)
    _tmp_cont = Continuum(test_ion, temperature_array)
    assert np.all(temperature_array == _tmp_cont.Temperature)


def test_free_free_scalar():
    # call free-free emission and loss rate methods for scalar temperature
    # emission--wavelength scalar
    tmp_cont_scalar.calculate_free_free_emission(wavelength_scalar)
    assert hasattr(tmp_cont_scalar, 'free_free_emission')
    assert tmp_cont_scalar.free_free_emission.shape == (1,1)
    # emission--wavelength array
    tmp_cont_scalar.calculate_free_free_emission(wavelength_array)
    assert hasattr(tmp_cont_scalar, 'free_free_emission')
    assert tmp_cont_scalar.free_free_emission.shape == (1,)+wavelength_array.shape
    # loss
    tmp_cont_scalar.calculate_free_free_loss()
    assert hasattr(tmp_cont_scalar, 'free_free_loss')
    assert tmp_cont_scalar.free_free_loss.shape == (1,)


def test_free_free_array():
    # call free-free emission and loss rate methods for scalar temperature
    # emission--wavelength scalar
    tmp_cont_array.calculate_free_free_emission(wavelength_scalar)
    assert hasattr(tmp_cont_array, 'free_free_emission')
    assert tmp_cont_array.free_free_emission.shape == temperature_array.shape+(1,)
    # emission--wavelength array
    tmp_cont_array.calculate_free_free_emission(wavelength_array)
    assert hasattr(tmp_cont_array, 'free_free_emission')
    assert tmp_cont_array.free_free_emission.shape == temperature_array.shape+wavelength_array.shape
    # loss
    tmp_cont_array.calculate_free_free_loss()
    assert hasattr(tmp_cont_array, 'free_free_loss')
    assert tmp_cont_array.free_free_loss.shape == temperature_array.shape


def test_free_bound_scalar():
    # free-bound emission and loss rate methods for scalar temperature
    # emiss--wavelength scalar
    tmp_cont_scalar.freeBound(wavelength_scalar)
    assert hasattr(tmp_cont_scalar, 'FreeBound')
    assert tmp_cont_scalar.FreeBound['intensity'].shape == (1,1)
    # emiss--wavelength array
    tmp_cont_scalar.freeBound(wavelength_array)
    assert hasattr(tmp_cont_scalar,'FreeBound')
    assert tmp_cont_scalar.FreeBound['intensity'].shape == (1,)+wavelength_array.shape
    # loss
    tmp_cont_scalar.calculate_free_bound_loss()
    assert hasattr(tmp_cont_scalar, 'free_bound_loss')
    assert tmp_cont_scalar.free_bound_loss.shape == (1,)


def test_free_bound_array():
    # free-bound emission and loss rate methods for temperature array
    # emiss--wavelength scalar
    tmp_cont_array.freeBound(wavelength_scalar)
    assert hasattr(tmp_cont_array, 'FreeBound')
    assert tmp_cont_array.FreeBound['intensity'].shape == temperature_array.shape+(1,)
    # emiss--wavelength array
    tmp_cont_array.freeBound(wavelength_array)
    assert hasattr(tmp_cont_array,'free_bound_emission')
    assert tmp_cont_array.FreeBound['intensity'].shape == temperature_array.shape+wavelength_array.shape
    # loss
    tmp_cont_array.calculate_free_bound_loss()
    assert hasattr(tmp_cont_array, 'free_bound_loss')
    assert tmp_cont_array.free_bound_loss.shape == temperature_array.shape


def test_free_bound_no_info():
    # raise error if no free-bound information is available
    with pytest.raises(ValueError, message='Expecting ValueError when no free-bound information is available'):
        tmp_cont_no_fb.freeBound(wavelength_array)
