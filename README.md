# ChiantiPy - Version 0.15.3
[![Documentation Status](http://readthedocs.org/projects/chiantipy/badge/?version=latest)](http://chiantipy.readthedocs.io/en/latest/?badge=latest)
[![Coverage Status](https://coveralls.io/repos/github/chianti-atomic/ChiantiPy/badge.svg?branch=master)](https://coveralls.io/github/chianti-atomic/ChiantiPy?branch=master)
[![ascl:1308.017](https://img.shields.io/badge/ascl-1308.017-blue.svg?colorB=262255)](http://ascl.net/1308.017)

ChiantiPy is the Python interface to the [CHIANTI atomic database](http://www.chiantidatabase.org) for astrophysical spectroscopy.  It provides the capability to calculate the emission line and continuum spectrum of an optically thin plasma based on the data in the CHIANTI database.

## What is CHIANTI?
CHIANTI provides a database of atomic data that can be used to interpret the emission of spectral lines and continua emitted from high-temperature, optically-thin astrophysical sources.  The CHIANTI project provides a suite of routines written in Interactive Data Language (IDL) to access the database and calculate various quantities for use in interpreting observed spectra or producing synthetic spectra.  As of ChiantiPy 0.10.0, the CHIANTI database version 10 or later is required

## Installation

ChiantiPy is available on [PyPI](https://pypi.org/project/ChiantiPy/)

```Shell
$ python3 -m pip install ChiantiPy
```

Or, clone and install the source from GitHub,
```Shell
$ git clone --recursive https://github.com/chianti-atomic/ChiantiPy.git
$ cd ChiantiPy
$ python setup.py install
```

Another method is to use the [Anaconda platform](https://www.continuum.io/downloads). Next, download the [CHIANTI database](http://www.chiantidatabase.org/chianti_download.html), version 10.0 or later. Assuming you've placed the CHIANTI tree in `$HOME`, set the environment variable in your `.bashrc` file,

```Shell
setenv XUVTOP $HOME/MY_CHIANTI_DIRECTORY

or

export XUVTOP=$HOME/MY_CHIANTI_DIRECTORY
```
should point to the top directory of your CHIANTI distribution



## Usage
As a quick example, we'll calculate the populations of the top 10 levels of Fe XIV as a function of temperature at constant density and plot them:
```Python
>>> import ChiantiPy.core as ch
>>> import numpy as np
>>> import matplotlib.pyplot as plt
>>> temperature = np.logspace(5.8,6.8,21)
>>> fe14 = ch.ion('fe_14',temperature=temperature,eDensity=1.e+9,em=1.e+27)
>>> fe14.popPlot()
>>> plt.show()
```

## Help
For more information about installing and using either ChiantiPy or the CHIANTI atomic database, check out the following links:

* [ChiantiPy Documentation on ReadTheDocs](https://chiantipy.readthedocs.io/)
* [ChiantiPy Documentation on ChiantiPy.github.io](https://chianti-atomic.github.io/)
* [Chianti Google Mailing List](https://groups.google.com/forum/#!forum/chianti)
* [CHIANTI Atomic Database Webpage](http://www.chiantidatabase.org/)
