<h1>ChiantiPy - Version 0.12.0</h1>
<p><a href="http://chiantipy.readthedocs.io/en/latest/?badge=latest"><img alt="Documentation Status" src="http://readthedocs.org/projects/chiantipy/badge/?version=latest" /></a>
<a href="https://coveralls.io/github/chianti-atomic/ChiantiPy?branch=master"><img alt="Coverage Status" src="https://coveralls.io/repos/github/chianti-atomic/ChiantiPy/badge.svg?branch=master" /></a>
<a href="http://ascl.net/1308.017"><img alt="ascl:1308.017" src="https://img.shields.io/badge/ascl-1308.017-blue.svg?colorB=262255" /></a></p>
<p>ChiantiPy is the Python interface to the <a href="http://www.chiantidatabase.org">CHIANTI atomic database</a> for astrophysical spectroscopy.  It provides the capability to calculate the emission line and continuum spectrum of an optically thin plasma based on the data in the CHIANTI database.</p>
<h2>What is CHIANTI?</h2>
<p>CHIANTI provides a database of atomic data that can be used to interpret the emission of spectral lines and continua emitted from high-temperature, optically-thin astrophysical sources.  The CHIANTI project provides a suite of routines written in Interactive Data Language (IDL) to access the database and calculate various quantities for use in interpreting observed spectra or producing synthetic spectra.  As of ChiantiPy 0.10.0, the CHIANTI database version 10 or later is required</p>
<h2>Installation</h2>
<p>The following dependencies are required to run ChiantiPy,</p>
<ul>
<li><a href="https://www.python.org/">Python3</a> (Python 3 is now required as of version 0.8.0)</li>
<li><a href="http://www.numpy.org/">Numpy</a></li>
<li><a href="https://www.scipy.org/">Scipy</a></li>
<li><a href="http://matplotlib.org/">Matplotlib</a></li>
<li><a href="https://github.com/ipython/ipyparallel">ipyparallel</a></li>
</ul>
<p>The following two are extremely useful for running Python programs
* <a href="http://ipython.org">IPython</a>
* <a href="http://jupyter.org/">Jupyter</a></p>
<p>Optionally, if you'd like to use the GUI dialogs,</p>
<ul>
<li><a href="https://riverbankcomputing.com/software/pyqt/intro">PyQt5</a></li>
</ul>
<p>If you are not familiar with installing Python and the needed dependencies, we recommend the <a href="https://www.continuum.io/downloads">Anaconda platform</a>. Next, download the <a href="http://www.chiantidatabase.org/chianti_download.html">CHIANTI database</a>, version 10.0 or later. Assuming you've placed the CHIANTI tree in <code>$HOME</code>, set the environment variable in your <code>.bashrc</code> file,
<code>Shell
export XUVTOP=$HOME/chianti/dbase</code></p>
<p>Finally, clone and install the source from GitHub,
<code>Shell
$ git clone --recursive https://github.com/chianti-atomic/ChiantiPy.git
$ cd ChiantiPy
$ python setup.py install</code>
The release is also available on <a href="https://pypi.org/project/ChiantiPy/">PyPI</a></p>
<h2>Usage</h2>
<p>As a quick example, we'll calculate the populations of the top 10 levels of Fe XIV as a function of temperature at constant density and plot them:
```Python</p>
<blockquote>
<blockquote>
<blockquote>
<p>import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt
temperature = np.logspace(5.8,6.8,21)
fe14 = ch.ion('fe_14',temperature=temperature,eDensity=1.e+9,em=1.e+27)
fe14.popPlot()
plt.show()
```</p>
</blockquote>
</blockquote>
</blockquote>
<h2>Help</h2>
<p>For more information about installing and using either ChiantiPy or the CHIANTI atomic database, check out the following links:</p>
<ul>
<li><a href="https://chiantipy.readthedocs.io/">ChiantiPy Documentation on ReadTheDocs</a></li>
<li><a href="https://groups.google.com/forum/#!forum/chianti">Chianti Google Mailing List</a></li>
<li><a href="http://www.chiantidatabase.org/">CHIANTI Atomic Database Webpage</a></li>
</ul>