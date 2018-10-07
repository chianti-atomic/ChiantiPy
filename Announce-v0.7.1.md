ChiantiPy version 0.7.1
=======================

changes from 0.7.0 to 0.7.1
---------------------------

version 0.7.0 included some changes in the ChiantiPy naming conventions, largely in the continuum class. These are being reverted to the original ChiantiPy naming conventions.

the ion.freeBoundxxx methods have been fixed and this also fixes the problem with the RadLoss class.

a pseudo-voigt filter has been added to tools.filters

the keyword argument wvlRange has been removed from the ion.emiss and ion.intensity methods

the keyword argument for the Emission Measure, em, has been removed from the ion.intensity and similar methods. It is now necessary to specify the emission when the object is instantiated.

a set of PyQt5 dialogs have been developed by **ktritz** and are now included

this is the last release that will use the PyQt4 widgets as an option.

the method **ioneqOne** is used by both the Ion and Continuum class. It has been moved to a single \_IoneqOne.py file in the **base** directory

current development is with Python 3.6

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
