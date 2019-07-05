===========
Changelog
===========

Changes from 0.8.6 to 0.8.7
===========================

continued code cleanup


Changes from 0.8.5 to 0.8.6
===========================

another bug-fix release

added argCheck method to make sure that sizes of temperature, density and emission measure were compatible

Changes from 0.8.4 to 0.8.5
===========================

This is a major bug-fix release.
================================

Errors in calculating the proton rates were corrected.

All temperatures and densities are now numpy arrays


Changes from 0.8.3 to 0.8.4
===========================

This is a major bug-fix release.
================================

Another significant bug was fixed in the important ion.populate method.


Changes from 0.7.1 to 0.8.3
===========================

This is a major bug-fix release.
================================

a small but mighty bug was found in the important ion.populate method.

Version 0.8.x files are necessary to use with the new CHIANTI Version 9.0 database
==================================================================================

Changes have been made to take into account the new way that CHIANTI is handling dielectronic recombination and autoionization

The release is also available on [PyPI](https://pypi.org/project/ChiantiPy/)

Documentation is available on [github.io](https://chianti-atomic.github.io/)

and on [ReadTheDocs](https://chiantipy.readthedocs.io/en/latest/?badge=latest)


changes from 0.7.1 to 0.8.0
===========================

ChiantiPy is now only compliant with Python 3.  Development is currently with Python 3.6

The use of the PyQt4 and WxWidgets packages have been dropped and PyQt5 is now used

The documentation is now available on github.io_ and ReadTheDocs_

.. _github.io:  https://chianti-atomic.github.io/

.. _ReadTheDocs:  https://chiantipy.readthedocs.io/en/latest/?badge=latest

changes from 0.7.0 to 0.7.1
===========================

version 0.7.0 included some changes in the ChiantiPy naming conventions, largely in the continuum class.  These are being reverted to the original ChiantiPy naming conventions.

the ion.freeBoundxxx methods have been fixed and this also fixes the problem with the RadLoss class.

a pseudo-voigt filter has been added to tools.filters

the keyword argument wvlRange has been removed from the ion.emiss and ion.intensity methods

the keyword argument for the Emission Measure, em, has been removed from the ion.intensity and similar methods.  It is now necessary to specify
the emission when the object is instantiated.

a set of PyQt5 dialogs have been developed by **ktritz** and are now included

this is the last release that will use the PyQt4 widgets as an option.

the method **ioneqOne** is used by both the Ion and Continuum class.  It has been moved to a single _IoneqOne.py file in the **base** directory


changes from 0.6.5 to 0.7.0
===========================

The primary change is that code developement has been moved to Github_.

.. _Github:  https://github.com/chianti-atomic/ChiantiPy

Also, in order to be more compliant with other astrophysical packages on Github (Astropy_ and SunPy_) the directory layout has been changed and renamed.


.. _Astropy:  https/github.com/astropy
.. _SunPy:  https://github.com/sunpy/sunpy

The core routines are now imported as

::

  import ChiantiPy.core as ch

this give access to ch.ion, ch.spectrum, etc.

In terms of bug-fixes, the calculation of excitation-autoionization cross-sections and rates have been corrected in the eaCross() and eaRate() methods

Current development is with Python 3.4

changes from 0.6.0 to 0.6.5
===========================

matplotlib.pyplot is now imported for plotting

IPython version 4 / Jupyter is now listed as a prerequisite.  However, v0.6.4 can be made compatible with IPython 2 or 3 with a simple edit.

An error in calculating the proton excitation rates was fixed.

The code has been edited to make it compatible with Python 3 and has been tested against Python 3.3

changes from 0.5.3 to 0.6.0
===========================

This is a major release.

First, ChiantiPy 0.6.0 is compatible with the most recently released CHIANTI database version 8.0.  It also fixes some major bugs in the previous version.  Documentation has been improved and a IPython notebook **QuickStart.ipynb**, that largely follows the 'Quick Start' documentation pages, has also been included.

There are two new multi-ion classes:  **bunch** and **ipymspectrum**.  **bunch** allows the user to calculate line intensities for a specified set of elements or individual ions as a function of temperature or density.  One advantage of **bunch** is the ability to calculate the intensity ratio of lines of two different ions as a function of temperature or density.

**ipymspectrum** is much like the existing **spectrum** and **mspectrum** classes.  **mspectrum** allows the use of the Python **multiprocessing** module to speed up spectral calculations.  The **ipymspectrum** class uses the IPython **parallel** module so that multiprocessing spectral calculations can be performed in the IPython QtConsole and Notebook.

A new method **intensityList** has been developed to allow the user to list the most intense lines within a given wavelength range.  This new methods, together with previously existing **intensityRatio** and **intensityRatioSave** are all now inherited by the **ion** classs and the  multi-ion classes.

The **ion** and multi-ion classes now accept the keyword argument **abundanceName** that allow the user to specify the set of elemental abundances rather than just the default abundance file.

Additional we have replaced the FortranFormat module of Scientific Python by Konrad Hinsen with the **fortranformat** module of Brendan Arnold at http://bitbucket.org/brendanarnold/py-fortranformat.  I have slightly modified fortranformat to make it Python 3 compliant.

For the future, I plan to make ChiantiPy compliant with both Python 2.7 and the current version of Python 3 (now 3.4), improve the documentation and move the projec to github, in no particular order.

ChiantiPy is now released under a new license, the OSI approved ISCL license.  From Wikipedia_ *The ISCL license is a permissive free software license written by the Internet Software Consortium (ISC). It is functionally equivalent to the simplified BSD and MIT/Expat licenses, ...*

.. _Wikipedia: https://en.wikipedia.org/w/index.php?title=ISC_license&oldid=664696993
