===================================
 Welcome to ChiantiPy version 0.6.4
===================================

changes from 0.6.0 to 0.6.4
===========================

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



What is ChiantiPy
=================

ChiantiPy is the Python interface to the CHIANTI atomic database for astrophysical spectroscopy.  It provides the capability to calculate the emission line and continuum spectrum of an optically thin plasma based on the data in the CHIANTI database.

Detailed information can be found at http://chiantipy.sourceforge.net

What is CHIANTI
===============

CHIANTI provides a database of atomic data that can be used to interpret the emission of spectral lines and continuua emitted from high-temperature, optically-thin astrophysical sources.  The CHIANTI project provides a suite of routines written in Interactive Data Language (IDL) to access the database and calculate various quantities for use in interpreting observed spectra or producing synthetic spectra.

==============================                                                                                                                
Getting started with ChiantiPy                                                                                                                
==============================                                                                                                                

Prerequisites
=============

* Python 2.7 or 3

* Numpy

* Scipy

* Matplotlib

* [Optional] PyQt4 or wxPython

*  IPython_ version 4 / Jupyter

  ChiantiPy has been developed with IPython versions 2.x, 3.x and now version 4 / Jupyter which is the required version.  A very small edit to one of the version 0.6.3 files will allow this version to work in IPython 2 or 3.    This is discussed in the Notes section.  The previously released version 0.6.0 is compatible with IPython 2.x and 3.x and will remain on the SourceForge site for an indefinite period of time.
  
.. _IPython:  http://ipython.org

* CHIANTI_, the atomic database for astrophysical spectroscopy

.. _CHIANTI: http://www.chiantidatabase.org

In addition, the FortranFormat module from Scientific Python, developed by Konrad Hinsen of the Centre de Biophysique Moleculaire (http://dirac.cnrs-orleans.fr/ScientificPython/), is included in this distribution for simplicity.

Installing the CHIANTI database
-------------------------------

The gzipped *data* tar ball can be downloaded from the CHIANTI website_

.. _website: http://www.chiantidatabase.org/download.html

*  put the file in a convenient directory, cd to the directory and untar the file

* ChiantiPy uses the environment variable *XUVTOP* to find the database.  Set XUVTOP to the name of the directory where the CHIANTI data tarball was placed.  For example

> setenv XUVTOP /data1/xuv/directory.where.the.tarball.was.placed

Some sites have the CHIANTI database maintained as part of a SolarSoft distribution.  In that case, simply set XUVTOP to the directory were it exists, usually something like $SSW/packages/chianti/dbase


Installing the ChiantiPy package
--------------------------------

Fairly detailed directions can be found on the web page ChiantiPy_

.. _ChiantiPy:  http://chiantipy.sourceforge.net/

The ChiantiPy package can be downloaded from the ChiantiPy, untar it, cd to the directory where it was unpacked, and then, as root

> python setup.py install

If you do not have root privileges, simply put the ChiantiPy directory (simply called 'chianti') in your PYTHONPATH

Version 0.6.4 is also available on the Python package index PyPI_

.. _PyPI: http://pypi.python.org


Running ChiantiPy
-----------------

The documentation can be found on its web page ChiantiPy_

.. _ChiantiPy:  http://chiantipy.sourceforge.net/

In particular, a quick start guide is included which should get you up and running fairly quickly.  This has been recently updated to reflect the new features in v0.6.0.


Keeping track of ChiantiPy
--------------------------

There is a mailing list that you can subscribe to at https://lists.sourceforge.net/lists/listinfo/chiantipy-users.  In order to subscribe it is first necessary to obtain a user account from sourceforge.net.  This is a straightforward process.

There is also a general chianti google group with the email address chianti@googlegroups.com
