import platform
import warnings

import numpy as np
import pytest

from ChiantiPy.core import spectrum, bunch, mspectrum, radLoss

try:
    from ChiantiPy.core import ipymspectrum
except ImportError as error:
    warnings.warn("Could not import IpyMspectrum; you may be missing ipyparellel.")
    ipymspectrum = error

# These tests verify that the spectrum classes (Spectrum, Mspectrum, and IpyMspectrum) and the Bunch class accept the
# same kind of basic array/container arguments for temperature, density, emission measure and minimum abundance.

_spectrum_like_type = (spectrum, bunch, mspectrum)
_container_type = (np.array,)  # tuple, list, np.atleast_2d)
_array_lengths = (1, 5)


# Test that scalar temperature and container works together.
@pytest.mark.parametrize("container_type", _container_type)
@pytest.mark.parametrize("density_length", _array_lengths)
def test_radloss_scalar_container_compatibility(container_type, density_length):
    temperature = 1e6
    density = container_type(np.linspace(1e9, 2e9, density_length))
    _tmp_spec = radLoss(temperature, density, minAbund=1e-3)


# Test that container and scalar density works together.
@pytest.mark.parametrize("container_type", _container_type)
@pytest.mark.parametrize("temperature_length", _array_lengths)
def test_radloss_container_scalar_compatibility(container_type, temperature_length):
    temperature = container_type(np.linspace(1e6, 2e6, temperature_length))
    density = 1e9
    _tmp_spec = radLoss(temperature, density, minAbund=1e-3)


# Test that the shape of the arguments are compatible. All arguments have to have length 1 on length N.
@pytest.mark.parametrize("container_type", _container_type)
@pytest.mark.parametrize("temperature_length", _array_lengths)
@pytest.mark.parametrize("density_length", _array_lengths)
def test_radloss_numpy_length_compatibility(container_type, temperature_length, density_length):
    temperature = container_type(np.linspace(1e6, 2e6, temperature_length))
    density = container_type(np.linspace(1e9, 2e9, density_length))
    _tmp_spec = radLoss(temperature, density, minAbund=1e-3)


# Save test as above but for the spectrum classes and bunch.
# The emission spectrum can also be None so that is also tested.
# Same test as above but also tests when all arguments are tuples or lists. Most likely these tests will succeed if
# the corresponding test_spectrum_numpy_length_compatibility succeeds.
@pytest.mark.parametrize("spectrum_like_type", _spectrum_like_type)
@pytest.mark.parametrize("container_type", _container_type)
@pytest.mark.parametrize("temperature_length", _array_lengths)
@pytest.mark.parametrize("density_length", _array_lengths)
@pytest.mark.parametrize("em_length", _array_lengths + (None,))
def test_spectrum_container_length_compatibility(spectrum_like_type, container_type,
                                                 temperature_length, density_length, em_length):
    if spectrum_like_type is mspectrum and platform.system() == 'Darwin':
        pytest.xfail("The %s class fails when running on MacOS." % str(spectrum_like_type))

    if type(spectrum_like_type) is ImportError:
        pytest.xfail("The %s class could not be imported; you may be missing ipyparallel." % str(spectrum_like_type))

    if em_length is not None:
        emission_measure = container_type(np.linspace(1e16, 2e16, temperature_length))
    else:
        emission_measure = None

    temperature = container_type(np.linspace(1e6, 2e6, temperature_length))
    density = container_type(np.linspace(1e9, 2e9, density_length))

    wavelength = np.linspace(200, 400)  # Just a few wavelengths.

    _tmp_spec = spectrum_like_type(temperature, density, wavelength, em=emission_measure, minAbund=1e-3)


# This tests that scalar values are also acceptable for density, temperature, and emission measure.
@pytest.mark.parametrize("spectrum_like_type", _spectrum_like_type)
@pytest.mark.parametrize("container_type", _container_type)
@pytest.mark.parametrize("temperature_length", _array_lengths)
@pytest.mark.parametrize("density_length", _array_lengths)
@pytest.mark.parametrize("em_length", _array_lengths + (None,))
def test_spectrum_container_scalar_compatibility(spectrum_like_type, container_type,
                                                 temperature_length, density_length, em_length):
    if spectrum_like_type is mspectrum and platform.system() == 'Darwin':
        pytest.xfail("The %s class fails when running on MacOS." % str(spectrum_like_type))

    if type(spectrum_like_type) is ImportError:
        pytest.xfail("The %s class could not be imported; you may be missing ipyparallel." % str(spectrum_like_type))

    if em_length is not None:
        emission_measure = container_type(np.linspace(1e16, 2e16, density_length))
    else:
        emission_measure = None

    temperature = container_type(np.linspace(1e6, 2e6, temperature_length))
    density = container_type(np.linspace(1e9, 2e9, density_length))

    wavelength = np.linspace(200, 400)  # Just a few wavelengths.

    # Temperature is scalar
    _tmp_spec = spectrum_like_type(1e6, density, wavelength, em=emission_measure, minAbund=1e-3)

    # Density is scalar
    _tmp_spec = spectrum_like_type(temperature, 1e9, wavelength, em=emission_measure, minAbund=1e-3)

    # Emission measure is scalar
    _tmp_spec = spectrum_like_type(temperature, density, wavelength, em=1e16, minAbund=1e-3)
