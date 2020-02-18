This directory contains 4 Jupyter IPython notebooks

QuickStart.ipyn - a notebook that generally follows the Quick-Start guide in the docs

The directory also contains 3 other notebooks and a json file.  These are demo files for reproducing some of the analyses in the paper "Electron densities and their uncertainties derived from spectral emission line intensities" by Kenneth Dere.  This paper has been submitted to the Monthly Notices of the Royal Astronomical Society.

There are two notebook files, fe_13_demo_chi2.ipynb and fe_13_demo_mcmc.ipynb for calculating the electron densities.  The first, by means of a Chi-squared minimization procedure and the second by means of Monte-Carlo-Markov-chain Bayesian inference using the PyMC package (https://github.com/pymc-devs/pymc).

The notebook file 'fe_13_demo_make_model.ipynb' constructs the model that is used by fe_13_demo_chi2.ipynb and fe_13_demo_mcmc.ipynb by reading the 'tab2_1993_qs_fe_13.json' file.  This files contains the Fe XIII line intensities from Brosius et al., 1996, Astrophysical Journal Supplement Series, 106, 143.  The notebook file needs to be run first.


