'''
ChiantiPy - CHIANTI Python package
Calculates various aspects of emission line and continua from the CHIANTI atomic database for astrophysical spectroscopy.
This is not yet an Astropy affiliated package, but it makes use of the Astropy package template
'''
# Include this stuff for astropy affiliation
# -------------------------------------------
from ._astropy_init import *
#--------------------------------------------

#if not _ASTROPY_SETUP_:
#    # For ChiantiPy
#    from . import  version
#    Version = version._last_generated_version
#
try:
    from .version import version as __version__
except ImportError:
    __version__ = '0.7.0'
try:
    from .version import githash as __githash__
except ImportError:
    __githash__ = ''
