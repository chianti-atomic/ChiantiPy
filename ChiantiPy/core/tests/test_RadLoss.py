"""
Tests for the radloss class
"""

import pytest
import numpy as np

from ChiantiPy.core import radLoss

temperature = np.logspace(5,8,20)
density = 1e+9
minAbund = 1e-4


# FIXME: should probably look at more cases
def test_radloss():
    tmp_radloss = radLoss(temperature, density, minAbund=minAbund)
    assert all((k in tmp_radloss.RadLoss for k in ('rate', 'temperature', 'density', 'minAbund', 'abundance')))

# TODO: test plotting
