"""
Tests for the ion class.
"""

import numpy as np
import pytest

from ChiantiPy.core import ion
import ChiantiPy.tools as ch_tools

# use an ion with relatively small chianti files
test_ion = 'fe_15'
# TODO: probably should check a few different ions, i.e. dielectronic, some without certain kinds
# files, etc.
# set a few temperature and density arrays
temperature_1 = 1e+6
temperature_2 = np.logspace(5,8,20)
density_1 = 1e+9
density_2 = np.logspace(5,8,20)
density_3 = np.logspace(5,8,21)
# setup an ion object to reuse in several tests
tmp_ion = ion(test_ion, temperature=temperature_2, eDensity=density_2)


# Check various ways to specify the temperature and density
def test_temperature_density():
    # TODO: test case where neither are set/ one or the other is not set
    # Two single values
    _tmp_ion = ion(test_ion, temperature=temperature_1, eDensity=density_1, setup=True)
    assert _tmp_ion.Temperature == np.array(temperature_1)
    assert _tmp_ion.EDensity == np.array(density_1)
    # Multiple temperatures, one density
    _tmp_ion = ion(test_ion, temperature=temperature_2, eDensity=density_1, setup=True)
    assert np.all(_tmp_ion.Temperature == temperature_2)
    assert np.all(_tmp_ion.EDensity == np.array(temperature_2.size*[density_1]))
    # One temperature, multiple densities
    _tmp_ion = ion(test_ion, temperature=temperature_1, eDensity=density_2, setup=True)
    assert np.all(_tmp_ion.Temperature == np.array(density_2.size*[temperature_1]))
    assert np.all(_tmp_ion.EDensity == density_2)
    # Two equal-sized temperature and density arrays
    _tmp_ion = ion(test_ion, temperature=temperature_2, eDensity=density_2, setup=True)
    assert np.all(_tmp_ion.Temperature == temperature_2)
    assert np.all(_tmp_ion.EDensity == density_2)
    # Two unequal sized temperature and density arrays
    with pytest.raises(ValueError,
                       message='''Expecting ValueError when temperature and density are not of
                                equal size.'''):
        _tmp_ion = ion(test_ion, temperature=temperature_2, eDensity=density_3, setup=True)


# Check how abundance is set
def test_abundance():
    # Float value
    _tmp_ion = ion(test_ion, temperature = temperature_1, eDensity = density_1, abundance=0.01, setup=False)
    assert _tmp_ion.Abundance == 0.01
    # FIXME: if setting custom abundance, AbundanceName should not be set, but right
    # now it is by the proton/electron density ratio calculation.
    # Custom filename
    _tmp_ion = ion(test_ion, temperature = temperature_1, eDensity = density_1, abundance='sun_coronal_2012_schmelz', setup=False)
    assert _tmp_ion.AbundanceName == 'sun_coronal_2012_schmelz'
    abundance = ch_tools.io.abundanceRead(abundancename='sun_coronal_2012_schmelz')
    assert _tmp_ion.Abundance == abundance['abundance'][_tmp_ion.Z-1]


# Check CHIANTI file imports
def test_chianti_files():
    assert hasattr(tmp_ion, 'Elvlc')
    assert hasattr(tmp_ion, 'Wgfa')
    if tmp_ion.Nscups > 0:
        assert hasattr(tmp_ion, 'Scups')
    if tmp_ion.Ncilvl > 0:
        assert hasattr(tmp_ion, 'Cilvl')
    if tmp_ion.Nreclvl > 0:
        assert hasattr(tmp_ion, 'Reclvl')
    if tmp_ion.Npsplups > 0:
        assert hasattr(tmp_ion, 'Psplups')
