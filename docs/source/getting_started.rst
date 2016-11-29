==============================
Getting started with ChiantiPy
==============================

Prerequisites
-------------

* CHIANTI_, the atomic database for astrophysical spectroscopy (Version 8)

.. _CHIANTI: http://www.chiantidatabase.org/

* Python_ (developed with versions 2.7 and 3.3, 3.4)

.. _PYTHON:  http://www.python.org

* Numpy_ (developed with 1.10)

.. _Numpy:  http://www.scipy.org/

* Scipy_ (developed with 0.16)

.. _Scipy:  http://www.scipy.org/

* Matplotlib_ (developed with 1.5)

.. _Matplotlib:  http://matplotlib.sourceforge.net/

* Matplotlib requires a GUI library

  PyQt4_ or wxPython_ (not compatible with Python3) or PyGTK_
  
  Once one of these is installed, it must be set as the backend in your matplotlibrc file, e.g., backend:  Qt4Agg

.. _PyQt4: http://www.riverbankcomputing.co.uk/

.. _wxPython:  http://www.wxpython.org/

.. _PyGTK:  http://www.pygtk.org/

* IPython_ version 4 or 5 / Jupyter
  
.. _IPython:  http://ipython.org



Install the CHIANTI database
----------------------------

The gzipped *data* tar ball can be downloaded from the CHIANTI website_

.. _website: http://www.chiantidatabase.org/chianti_download.html

*  put the file in a convenient directory, cd to the directory and untar the file

* ChiantiPy uses the environment variable *XUVTOP* to find the database.  Set XUVTOP to the name of the directory where the CHIANTI data tarball was placed.  For example

::
	
  setenv XUVTOP /data1/xuv/directory.where.the.tarball.was.placed
  

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

In order to be compatible with the latest version (8.0) of the CHIANTI atomic database, it is necessary to install the latest version (0.6.0) of ChiantiPy

::

  pip install ChiantiPy


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

Note - ChiantiPy interactions with Matplotlib
---------------------------------------------------

Some of the ChiantiPy methods ask the user to make a selection.  With ChiantiPy, this can be done within the command line shell or with gui dialog widgets using PyQt4 or wxPython.  Matplotlib needs to have a backend specified.  The default is 'GTK' and ChiantiPy will use the command line shell for user input.  If the Matplotlib backend is specified to be 'Qt4Agg', then the PyQt4 widget set will be used by ChiantiPy.  If the Matplotlib backend is specified to be 'WXAgg' then the wxPython widget set will be used by ChiantiPy.  If the Matplotlib backend is set to something other than the 3 values previously discussed, the command line shell will be used.

