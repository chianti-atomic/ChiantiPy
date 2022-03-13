This directory contains 4 Jupyter IPython notebooks

QuickStart.ipynb - a notebook that generally follows the Quick-Start guide in the docs

The directory also contains 3 other notebooks and a json file.  These are demo files for reproducing some of the analyses in the paper "Electron densities and their uncertainties derived from spectral emission line intensities" by Kenneth Dere.  This paper has been published in the Monthly Notices of the Royal Astronomical Society, 2020, 496, 2334.

The notebook file '1_fe_13_demo_make_model.ipynb' constructs the model that is used by fe_13_demo_chi2.ipynb and fe_13_demo_mcmc.ipynb by reading the 'tab2_1993_qs_fe_13.json' file.  This files contains the Fe XIII line intensities from Brosius et al., 1996, Astrophysical Journal Supplement Series, 106, 143.  This notebook file needs to be run first.

The next notebook to run is:  2_fe_13_demo_check_model.ipynb  -- This notebook load the previously created pickle containing the match attribute.  In this notebook a density and emission measure are guessed from an 'em loci' plot and the predictions compard with the observations.

The next notebook to run is:  3_fe_13_demo_chi2_search.ipynb -- This performs a brute force chi-squared search over the density range and finds the best fit density and emission measure.  These best values are then inserted into the model, a prediction made and compared with the observations.

An other notebook files require revision.  This should happen soon.  This notebook is for determining the electron densities by means of Monte-Carlo-Markov-chain Bayesian inference using the PyMC3 package (https://github.com/pymc-devs/pymc).




