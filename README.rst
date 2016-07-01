===========
ChiantiPy
===========
.. image:: http://readthedocs.org/projects/chiantipy/badge/?version=latest
   :target: http://chiantipy.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://kiwiirc.com/buttons/chat.freenode.com/openastronomy.png
   :target: https://kiwiirc.com/client/chat.freenode.com/?nick=OpenAstron|?#openastronomy
   :alt: Visit the OpenAstronomy IRC channel

ChiantiPy is the Python interface to the CHIANTI atomic database for astrophysical spectroscopy.  It provides the capability to calculate the emission line and continuum spectrum of an optically thin plasma based on the data in the CHIANTI database.

CHIANTI_ provides a database of atomic data that can be used to interpret the emission of spectral lines and continuua emitted from high-temperature, optically-thin astrophysical sources.  The CHIANTI_ project provides a suite of routines written in Interactive Data Language (IDL) to access the database and calculate various quantities for use in interpreting observed spectra or producing synthetic spectra.

Installation
==============================
Dependencies
------------

* Python_ (2.7,3)
* Astropy_
* Numpy_
* Scipy_
* Matplotlib_
* IPython_ (>=v4)
* ipyparallel_

If you want to use the GUI:

* PyQt4_ and/or
* wxPython_

You will also need to install CHIANTI_, the atomic database for astrophysical spectroscopy. See below_ for more instructions on how to install the database.

In addition, the FortranFormat module from Scientific Python, developed by Konrad Hinsen of the `Centre de Biophysique Moleculaire <http://dirac.cnrs-orleans.fr/ScientificPython/>`_ is included in this distribution for simplicity.

.. _CHIANTI: http://www.chiantidatabase.org
.. _IPython:  http://ipython.org
.. _Python: https://www.python.org/
.. _Astropy: http://www.astropy.org/
.. _Numpy: http://www.numpy.org/
.. _Scipy: https://www.scipy.org/
.. _Matplotlib: http://matplotlib.org/
.. _ipyparallel: https://github.com/ipython/ipyparallel
.. _PyQt4: https://riverbankcomputing.com/software/pyqt/intro
.. _wxPython: http://www.wxpython.org/

Installing the CHIANTI database
-------------------------------
.. _below:

CHIANTI can be downloaded as a `gzipped data tar ball <http://www.chiantidatabase.org/download/CHIANTI_8.0.2_data.tar.gz>`_ or through `SSW <http://www.lmsal.com/solarsoft/sswdoc/sswdoc_jtop.html>`_. Detailed download instructions can be found `here <http://www.chiantidatabase.org/download.html>`_.

ChiantiPy uses the environment variable ``XUVTOP`` to find the database.  Set ``XUVTOP`` to the name of the directory where the CHIANTI database is located. For example, in bash,:

  export XUVTOP=/path/to/chianti/tree

Similarly, in csh,:

  setenv XUVTOP /path/to/chianti/tree

If you have installed CHIANTI via SSW, the database is probably located at ``/path/to/ssw/packages/chianti/dbase``.

Installing ChiantiPy
----------------------
The best way to obtain the ChiantiPy code is via GitHub. To clone this repository,:

  git clone --recursive https://github.com/chianti-atomic/ChiantiPy.git

The code in the ChiantiPy repository will change periodically as features are added and bugs are fixed. To pull down the latest changes, inside the directory where you cloned ChiantiPy, run:

  git pull

Alternatively, if you'd prefer not to use git or have not installed it, simply download the `zipped repository <https://github.com/chianti-atomic/ChiantiPy/archive/master.zip>`_.

After obtaining the code, ``cd`` into the directory where you installed ChiantiPy and run:

  python setup.py install

If you do not have root privileges, you can do a local install (e.g. in ``/path/to/local/install``) by first adding ``/path/to/local/install`` to your ``PYTHONPATH`` and then running:

    python setup.py install --prefix=/path/to/local/install

To test your installation, in a Python or IPython shell, run:

  import ChiantiPy.core as ch

If you don't get any import errors, you've successfully installed ChiantiPy!

While ChiantiPy is available via `PyPi <http://pypi.python.org>`_, v0.6.4 is the most current version uploaded there. For now, the most reliable way to obtain ChiantiPy is via GitHub.

Usage
=======
The full documentation can be found on the `Read the Docs <http://chiantipy.readthedocs.io/en/latest/>`_ page. The `quickstart guide <http://chiantipy.readthedocs.io/en/latest/quick_start.html>`_ shows several examples of the capabilities of ChiantiPy.

Additionally, a `Jupyter <http://jupyter.org/>`_ notebook showing some of these features can be found in ``ipython_notebooks/QuickStart.ipynb``.

Keeping track of ChiantiPy
===========================
There is a mailing list that you can subscribe to at https://lists.sourceforge.net/lists/listinfo/chiantipy-users.  In order to subscribe it is first necessary to obtain a user account from sourceforge.net.  This is a straightforward process.

There is also a general chianti google group with the email address chianti@googlegroups.com
