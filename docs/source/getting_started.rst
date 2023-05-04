==============================
Getting started with ChiantiPy
==============================

Prerequisites
-------------

* CHIANTI_, the atomic database for astrophysical spectroscopy (Version 10 or later)

.. _CHIANTI: http://www.chiantidatabase.org/

* Python3_ (3.8 is the current development version)

.. _PYTHON3:  http://www.python.org

* Numpy_ (currently developed with 1.20)

.. _Numpy:  http://www.scipy.org/

* Scipy_ (currently developed with 1.6)

.. _Scipy:  http://www.scipy.org/

* Matplotlib_ (currently developed with 3.3.4)

.. _Matplotlib:  http://matplotlib.sourceforge.net/

* Matplotlib requires a GUI library

  PyQt5_

  Once one of these is installed, it must be set as the backend in your matplotlibrc file, e.g., backend:  Qt5Agg

.. _PyQt5: http://www.riverbankcomputing.co.uk/


* ipyparallel (required for multiprocessing with ipymspectrum)

* (not really a prerequisite but **extremely** useful) IPython_ version 7.21 and Jupyter_

.. _IPython:  http://ipython.org

.. _Jupyter: http://jupyter.readthedocs.io/en/latest/


Install the CHIANTI database
----------------------------

The gzipped *data* tar ball can be downloaded from the CHIANTI website_

.. _website: http://www.chiantidatabase.org/chianti_download.html

*  put the file in a convenient directory, cd to the directory and untar the file

* ChiantiPy uses the environment variable *XUVTOP* to find the database.  Set XUVTOP to the name of the directory where the CHIANTI data tarball was placed.  For example

::

  setenv XUVTOP /data1/directory.where.the.tarball.was.placed


or on Windows:   To set the environment variable, go to Control Panel -> System -> Advanced System Properties -> Environment Variables.


Some sites have the CHIANTI database maintained as part of a SolarSoft distribution.  In that case, simply set XUVTOP to the directory were it resides, usually something like $SSW/packages/chianti/dbase

Install the Prerequisites
-------------------------

On **Linux** systems this can usually be done with your package manager.

On **Windows**, **Linux** and **Mac** systems, it is possible to use

* the Anaconda_ distribution from from Continuum, or,

.. _Anaconda:  http://continuum.io/downloads

* the Canopy_ distribution from Enthought.

.. _Canopy:  https://store.enthought.com/downloads/#default

On **Windows**, it is also possible to use:

* WinPython_

.. _WinPython:  http://winpython.github.io/

All of these packages are free, at least for noncommercial use (I believe) and have a considerable amount of documentation.  **You shoud check the version of IPython that is provided**.


Install the ChiantiPy package
-----------------------------

In order to be compatible with the latest version (10) of the CHIANTI atomic database, it is necessary to install the latest version (0.10.0) of ChiantiPy

::

  pip install ChiantiPy

or

::

  pip3 install ChiantiPy


I have not tried this with ChiantiPy, myself.


The ChiantiPy package can be downloaded from the ChiantiPy_ project page at Sourceforge, untar it, cd to the directory where it was unpacked, and then, as root

.. _ChiantiPy:  http://sourceforge.net/projects/chiantipy/

::

  python setup.py install

If you do not have root privileges, simply put the ChiantiPy directory in your PYTHONPATH

::

  python setup.py install --prefix=somewhere_in_my_PYTHONPATH


or on a Mac, with the Anaconda package

::

  python setup.py install --prefix=/Users/your_user_name/anaconda/

**Thanks to Peter Young (GMU) for providing the instructions for installation on Mac and Windows**

Note - ChiantiPy interactions with various GUI backends
-------------------------------------------------------

First, Matplotlib requires a GUI backend and can be specified by the user in the matplotlibrc file.  Matplotlib expects to find this file in  $HOME/.config/matplotlib, although it might require that you copy it to that directory.

ChiantiPy also uses a GUI dialog widget set.  Selections can also be made via the command line.  The user choice is specified in the chiantirc file.  One is included with the distribution.  On Linux, if it is copied into either the $HOME/.config or $HOME/.chianti directory, it will automatically be picked up.  It is a text file and can be edited.  On Windows, it should be copied to $PROFILEHOME/.config or $PROFILEHOME/.chianti where it will also be picked up.  Otherwise, the default GUI is to use the command line.  A default chiantirc file is included with the distribution.

In order for the ChiantiPy dialog widget to be used, a backend for them must be initiated.  If you choose the same backend for matplotlib (PyQt5 is best) as for the ChiantiPy widgets, then running %matplotlib in an IPython or Jupyter session will do the trick.  In an interactive Python session, invoking a matplotlib or matplotlib qt command first should also do the trick.

::

  matplotlib inline

or

::

  matplotlib qt

if you are using Qt5

If you choose to use a GUI backend other than that used for matplotlib, then in an IPython or a Jupyter command the following magic commands are also available to start the backend:

::

  %gui qt


ChiantiPy has mostly been tested with the Qt5 backend for Matplotlib and using the ChiantiPy Qt5 widgets.

