ChiantiPy - Version 0.8.5
====================================================================================

| |Documentation Status|
| |Coverage Status|
| |ascl:1308.017|

ChiantiPy is the Python interface to the `CHIANTI atomic
database <http://www.chiantidatabase.org>`__ for astrophysical
spectroscopy. It provides the capability to calculate the emission line
and continuum spectrum of an optically thin plasma based on the data in
the CHIANTI database.

What is CHIANTI?
----------------

CHIANTI provides a database of atomic data that can be used to interpret
the emission of spectral lines and continua emitted from
high-temperature, optically-thin astrophysical sources. The CHIANTI
project provides a suite of routines written in Interactive Data
Language (IDL) to access the database and calculate various quantities
for use in interpreting observed spectra or producing synthetic spectra.

Installation
------------

The following dependencies are required to run ChiantiPy,

-  `Python3 <https://www.python.org/>`__ (Python 3 is now required as of
   version 0.8.0)
-  `Numpy <http://www.numpy.org/>`__
-  `Scipy <https://www.scipy.org/>`__
-  `Matplotlib <http://matplotlib.org/>`__
-  `ipyparallel <https://github.com/ipython/ipyparallel>`__

The following two are extremely useful for running Python programs

-  `IPython <http://ipython.org>`__
-  `Jupyter <http://jupyter.org/>`__

Optionally, if you'd like to use the GUI dialogs,

-  `PyQt5 <https://riverbankcomputing.com/software/pyqt/intro>`__

If you are not familiar with installing Python and the needed
dependencies, we recommend the `Anaconda
platform <https://www.continuum.io/downloads>`__. Next, download the
`CHIANTI
database <http://www.chiantidatabase.org/chianti_download.html>`__,
version 9.0 or later. Assuming you've placed the CHIANTI tree in
``$HOME``, set the environment variable in your ``.bashrc`` file,

.. code:: shell

    export XUVTOP=$HOME/chianti/dbase

Finally, clone and install the source from GitHub,

.. code:: shell

    $ git clone --recursive https://github.com/chianti-atomic/ChiantiPy.git
    $ cd ChiantiPy
    $ python setup.py install

The release is also available on
`PyPI <https://pypi.org/project/ChiantiPy/>`__

Usage
-----

As a quick example, we'll calculate the populations of the top 10 levels
of Fe XIV as a function of temperature at constant density and plot
them:

.. code:: python

    >>> import ChiantiPy.core as ch
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> temperature = np.logspace(5.8,6.8,21)
    >>> fe14 = ch.ion('fe_14',temperature=temperature,eDensity=1.e+9,em=1.e+27)
    >>> fe14.popPlot()
    >>> plt.show()

Help
----

For more information about installing and using either ChiantiPy or the
CHIANTI atomic database, check out the following links:

-  `ChiantiPy Documentation on
   github.io <http://chianti-atomic.github.io/>`__
-  `ChiantiPy Documentation on
   ReadTheDocs <https://chiantipy.readthedocs.io/>`__
-  `Chianti Google Mailing
   List <https://groups.google.com/forum/#!forum/chianti>`__
-  `CHIANTI Atomic Database Webpage <http://www.chiantidatabase.org/>`__

.. |Documentation Status| image:: http://readthedocs.org/projects/chiantipy/badge/?version=latest
   :target: http://chiantipy.readthedocs.io/en/latest/?badge=latest
.. |Coverage Status| image:: https://coveralls.io/repos/github/chianti-atomic/ChiantiPy/badge.svg?branch=master
   :target: https://coveralls.io/github/chianti-atomic/ChiantiPy?branch=master
.. |ascl:1308.017| image:: https://img.shields.io/badge/ascl-1308.017-blue.svg?colorB=262255
   :target: http://ascl.net/1308.017
