ChiantiPy version 0.7.0
=======================

changes from 0.6.5 to 0.7.0
---------------------------

The is the first release of ChiantiPy on [Github](https://github.com/chianti-atomic/ChiantiPy).

This release is also available on [PyPI](https://pypi.python.org/pypi).

To be more compliant with other astrophysical packages on Github ([Astropy](http://astropy.org/) and [SunPy](https://github.com/sunpy/sunpy)) the directory layout has been changed and renamed.

The core routines are now imported as

    import ChiantiPy.core as ch

this give access to ch.ion, ch.spectrum, etc.

Bug-fixes: the calculation of excitation-autoionization cross-sections and rates have been corrected in the eaCross() and eaRate() methods

Current development is with Python 3.4
