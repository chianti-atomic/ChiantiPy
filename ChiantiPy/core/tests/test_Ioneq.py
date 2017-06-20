"""
Tests for ioneq class
"""

import numpy as np
import pytest

import ChiantiPy.tools.util as util
from ChiantiPy.core import ioneq

temperature = np.logspace(4,9,50)
# Element and Z should represent the same element
el = 'Fe'
z = 26

# TODO: test temperatures out of valid extrapolation range


def test_el_input():
    el_ioneq_input = ioneq(el)
    assert el_ioneq_input.Z == util.el2z(el)


def test_z_input():
    z_ioneq_input = ioneq(z)


def test_el_z_inputs_same():
    el_ioneq_input = ioneq(el)
    z_ioneq_input = ioneq(z)
    el_ioneq_input.load()
    z_ioneq_input.load()
    assert np.all(el_ioneq_input.Temperature == z_ioneq_input.Temperature)
    assert np.all(el_ioneq_input.Ioneq == z_ioneq_input.Ioneq)


def test_calculate_ioneq():
    calc_ioneq = ioneq(el)
    calc_ioneq.calculate(temperature)
    assert hasattr(calc_ioneq, 'Temperature')
    assert hasattr(calc_ioneq, 'Ioneq')
    assert np.all(calc_ioneq.Temperature == temperature)


def test_load_ioneq():
    load_ioneq = ioneq(el)
    load_ioneq.load()
    assert hasattr(load_ioneq, 'Temperature')
    assert hasattr(load_ioneq, 'Ioneq')


def test_load_ioneq_alternate_file():
    load_ioneq = ioneq(el)
    load_ioneq.load(ioneqName='chianti')
    assert hasattr(load_ioneq, 'Temperature')
    assert hasattr(load_ioneq, 'Ioneq')
