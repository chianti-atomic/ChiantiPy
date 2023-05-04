===========
Changelog
===========

Changes from 0.15.0 to 0.15.1
=============================

The continuum.freeBound method has been corrected so that it can be used under extreme conditions and not fail.

The io.abundanceRead method has been updated to work with new files in the upcoming CHIANTI 10.1 release


Changes from 0.14.1 to 0.15.0
=============================

Significant improvement have been made.  It is now possible to calculate spectral line intensities in wavelength/energy units of angstroms, nm, eV, or keV

Calculations of the continuum still require that the wavelengths are in angstroms.

The free-bound continuum now includes the photoionization cross sections of Verner for recombination to the ground level

The free-free continuum has been correction to use the ion charge, not the nuclear charge previously used

A new class MradLoss.mradLoss has been created.  It allows multiprocessing calculations of the radiative loss rate

The chiantirc file can now also be placed in the $HOME/.config directory


Changes from 0.14.0 to 0.14.1
=============================

This relatively minor release adds some new features and corrects some glitches

A function demRead has been to ChiantiPy.tools.io for reading the CHIANTI .dem files in the XUVTOP/dem directory

The spectrumPlot method has been updated to provide more correct labeling of syntheic spectra

The QuickStart guide (html and notebook) have been updated to reflect these changes and show how to use the .dem files


Changes from 0.13.1 to 0.14.0
=============================

a new class 'ch.redux' restores the attributes saved by the saveData methods.  It inherits as number of methods, especially, for plotting.

the inherited methods 'intensityPlot' and 'spectrumPlot' have been improved.  These are inherited by the ion, bunch, spectrum, mspectrum, ipymspectrum and redux classes.

First, these methods will also display the ion name ('Fe XIV') and wavelength together with the line intensities or spectral intensity.

These are more flexible and several keyword arguments have been added:

    'doLabel' governs whether to display the ion name and wavelength
    'lw' the linewidth of the marker in matplotlib units (default=1)
    'doTitle' governs whether to add a title to the plot

The QuickStart (html and notebook) has been updated to demonstrate some of these new features

Import bug fix:  the indices for calculating the two-photon continue were update to match the new ordering of the energy levels for the h-like and he-like ions.

The ionization potential array (Ip) has been enlarged to make room for call to zn_31

Changes from 0.13.0 to 0.13.1
=============================

This is primarily a bug fix release to correct a bug in ionGate that was not taken care of in the 0.13.0 release.


Changes from 0.12.0 to 0.13.0
=============================

it is now possible to incorporate a user created abundance file.  It needs to be of the same structure as one of the .abund files in the XUVTOP/abundance directory.  The file can be located anywhere on the user's computer.  It will be read if the abundance keyword is set to a fully qualified file name of the new abundance file, such as '/home/me/myabundance.abund', or equally, '/home/me/myabundance.txt'

in the core classes, bunch and spectrum, now have saveData methods to save calculations to a pickle file

by default, the pyQt widgets are now used a selection tools.  The command line selection tools are still available but the ~/.chianti/chiantirc file needs to select them.

Problems with the ionGate method have been fixed.  This method help to select which ions will be used in multi-ion calculations




Changes from 0.11.0 to 0.12.0
=============================


The model module is more mature

For Windows users, it is now possible to place the chiantirc file in $PROJECTHOME/.config or $PROJECTHOME/.chianti

Many improved docstrings for the documentation

the bunch class has been moved to a new module core.Bunch

a number of jupyter ipython notebooks have been created/improved to demonstrate the use of the bunch, spectrum and model.Maker classes.  A short README.txt can be found in the same directory provides an introduction to these notebooks

A bug in the inherited method base._IntensityRatio() was not properly corrected in v0.10.0.  This is fixed here



Changes from 0.10.0 to 0.11.0
=============================

Calculations of the free-bound/radiative recombination continuum and radiation losses depend on a file that provides and LS description of the bound and singly excited energy levels.  This file is called c_5.fblvl in the case of C V.  Not all ions have an associated .fblvl files and is was necessary to revise ChiantiPy to ignore the free-bound calculation for these ions.

In addition, it was found that under extreme conditions, such as very low temperatures for highly ionized species, that bad values would arise (Nans and infinities).  These are now detected and removed.


Changes from 0.9.5 to 0.10.0
============================

The Karzas and Latter (1965) bound-free gaunt factors in the CHIANTI database have been corrected as of CHIANTI version 10.  This effects the calculation of the free-bound continuum.  The continuum.freeBound method has been updated to uses these new data.

The freeBound and the freeFree methods now have 2 new keyword variables:  **includeAbund** and **includeIoneq**.  Their initial values are True so that elemental abundance and the ionization equilibrium appropriate to the ion is included the the output spectrum.

The freeBoundEmiss is removed as it has become redundant

A new continuum method, freeBoundLossMao includes the radiative-recombination (free-bound) loss rate as calculated by Mao et al. (2017)

a bug in the inherited method base._IntensityRatio() was corrected.



Changes from 0.9.4 to 0.9.5
===========================

this is a bug-fix release.

a bug in the inherited method base._IntensityRatio() had a problem if lines were selected from different ions


Changes from 0.9.3 to 0.9.4
===========================

this is a bug-fix release.

changes in version 0.9.2 continued to give problems with ions that included autoionization rates


Changes from 0.9.3 to 0.9.3
===========================

this is a bug-fix release.

changes in version 0.9.2 gave problems with ions that included autoionization rates



Changes from 0.9.2 to 0.9.3
===========================

this is a bug-fix release.

changes in version 0.9.2 led to an error where ion.recombRate did not work.  This has been fixed


Changes from 0.9.1 to 0.9.2
===========================

this is a bug-fix release.

changes in version 0.9.1 lead to an error where a bare ion has not recombination rate.  This has been fixed


Changes from 0.9.0 to 0.9.1
===========================

this is a bug-fix release.

in cases when it is not possible to calculate the free-bound continuum for some ion, mspectrum did not handle this correctly and crashed

also, the ion zn_31 (Zn XXXI) is a bare ion and has no ionization potential (IP) and looking it up caused indexing errors.


Changes from 0.8.7 to 0.9.0
===========================

a new module model.maker has been added

::

  import ChiantiPy.model as mdl
  mymodel = mdl.maker(...)


a serious bug in ch.freeBound was fixed - the use of a single temperature was problematic

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

The primary change is that code development has been moved to Github_.

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

For the future, I plan to make ChiantiPy compliant with both Python 2.7 and the current version of Python 3 (now 3.4), improve the documentation and move the project to github, in no particular order.

ChiantiPy is now released under a new license, the OSI approved ISCL license.  From Wikipedia_ *The ISCL license is a permissive free software license written by the Internet Software Consortium (ISC). It is functionally equivalent to the simplified BSD and MIT/Expat licenses, ...*

.. _Wikipedia: https://en.wikipedia.org/w/index.php?title=ISC_license&oldid=664696993
