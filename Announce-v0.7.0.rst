=======================
ChiantiPy version 0.7.0
=======================


changes from 0.7.0 to 0.7.1
===========================

version 0.7.0 included some changes in the ChiantiPy naming conventions, largely in the continuum class.  These are being reverted to the original ChiantiPy naming conventions.

a pseudo-voigt filter has been added to tools.filters

the keyword argument wvlRange has been removed from the ion.emiss and ion.intensity methods

a set of PyQt5 dialogs have been developed by **ktritz** and are now included

this is the last release that will use the PyQt4 widgets as an option.

Current development is with Python 3.6


changes from 0.6.5 to 0.7.0
===========================

The is the first release of ChiantiPy on Github_.

.. _Github:  https://github.com/chianti-atomic/ChiantiPy

This release is also available on PyPI_.

.. _PyPI:  https://pypi.python.org/pypi

To be more compliant with other astrophysical packages on Github (Astropy_ and SunPy_) the directory layout has been changed and renamed.

.. _Astropy:  http://astropy.org/

.. _SunPy:  https://github.com/sunpy/sunpy

The core routines are now imported as

::

  import ChiantiPy.core as ch

this give access to ch.ion, ch.spectrum, etc.

Bug-fixes:  the calculation of excitation-autoionization cross-sections and rates have been corrected in the eaCross() and eaRate() methods

Current development is with Python 3.4
