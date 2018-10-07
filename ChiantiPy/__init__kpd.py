"""
ChiantiPy - CHIANTI Python package Calculates various aspects of emission lines
and continua from the CHIANTI atomic database for astrophysical spectroscopy.
"""
# This is not yet an Astropy affiliated package, but it makes use of the Astropy
# package template

# this indicates whether or not we are in the package's setup.py
try:
    _ASTROPY_SETUP_
except NameError:
    from sys import version_info
    if version_info[0] >= 3:
        import builtins
    else:
        import __builtin__ as builtins
    builtins._ASTROPY_SETUP_ = False

try:
    from .version import version as __version__
except ImportError:
    __version__ = ''
try:
    from .version import githash as __githash__
except ImportError:
    __githash__ = ''

# Import astropy test runner if we can and dummy it if we can't
import os
try:
    from astropy.tests.helper import TestRunner
    test = TestRunner.make_test_runner_in(os.path.dirname(__file__))
except ImportError:
    def test(*args, **kwargs):
        raise ImportError("astropy is needed to run the tests")

# Actual package imports here:
# Note this if statement is only here to allow chiantipy to be imported before
# it's compiled.
if not _ASTROPY_SETUP_:
    ## For ChiantiPy
    #from . import version
    #Version = version._last_generated_version
    ## For ChiantiPy
    from . import  version
    #Version = version._last_generated_version
    __version__ = version.__version__
    __version_info__ = version.__version_info__
