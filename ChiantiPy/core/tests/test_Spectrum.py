"""
Tests for the spectrum and bunch classes
"""

import numpy as np
import pytest

from ChiantiPy.core import spectrum, bunch

# set temperature, density, wavelength
temperature_scalar = 2e6
temperature_array_1d = np.array([temperature_scalar])
temperature_array = np.array([1e+6,4e+6,1e+7])
temperature_array_long = np.logspace(5,8,10)
density = 1e+9
em = 1.e+27
wavelength = np.linspace(200,400,1000)
wavelength_range = [wavelength[0], wavelength[-1]]
min_abund = 1.e-4
ion_list = ['fe_15', 'fe_16']


def test_spectrum_scalar():
    _tmp_spec = spectrum(temperature_scalar, density, wavelength, minAbund=min_abund)
    assert _tmp_spec.Spectrum['intensity'].shape == wavelength.shape


def test_spectrum_array():
    _tmp_spec = spectrum(temperature_array_1d, density, wavelength, minAbund=min_abund)
    assert _tmp_spec.Spectrum['intensity'].shape == wavelength.shape
    _tmp_spec = spectrum(temperature_array, density, wavelength, minAbund=min_abund)
    assert _tmp_spec.Spectrum['intensity'].shape == temperature_array.shape+wavelength.shape


def test_bunch():
    _tmp_bunch = bunch(temperature_array_long, density, wvlRange=wavelength_range, ionList=ion_list)
    # TODO: need to assert something here, not clear what exactly to test yet
