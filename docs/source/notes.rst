=====
Notes
=====

Setting default values
----------------------

Several parameter values can be set to define the way ChiantiPy behaves.  The following parameter values can be set in your chiantirc file.  A chiantirc file is available in the downloaded package but it will only be read if it is located in either ~HOME/.config (probably the best place) or ~HOME/.chianti or ~PROFILEHOME/.config or ~PROFILEHOME/.chianti.  It is then possible to modifiy the file to select the default value you wish to use:

These are:


wavelength

    The default value is angstrom.  Other possible values are *nm* for nanometers, *ev* for electron volts or *kev* for kilo-electron volts.  However, the continuum classes can only be used with the wavelengths in angstroms


flux
    The default value is energy. Acceptable values are *energy* and *photon* and these govern emissivities and intensities.  If set to energy, emission units are in erg.

abundfile

    The name of the abundance file.  Acceptable values are any of the file names in XUVTOP/abundance, such as *cosmic_1973_allen*.  The default value is *sun_photospheric_2015_scott* that includes the abundances of Scott et al., 2015, A&A, 573, A25.

ioneqfile
    the name of the ionization equilibrium file.  Acceptable values are any of the file names in XUVTOP/ioneq such as *arnaud_raymond*, *arnaud_rothenflug*, or *chianti*.  The default value is *chianti* which includes the ionization equilibrium calculations of K.P. Dere, G. Del Zanna, P.R. Young, E. Landi, R Sutherland, 2019, ApJ, 241, 2 and are considered to be based on the best ionization and recombination rates currently available.



Setting *minAbund* in spectrum calculations
-------------------------------------------

When calculation spectra with *spectrum* or *mspectrum*, it is often useful to set the "minAbund" keyword which governs the minimum abundance of any element included in the calculation.  Below is a list of elemental abundances for the elements through zinc and the elements that will be included by several value of "minAbund".  These are for photospheric abundances, some of which may have changed since this table was made.

=======  =========  ========  =======  ======


Element  Abundance          minAbund
-------  ---------  -------------------------
  ..       ..       1.e-6     1.e-5     1.e-4
=======  =========  ========  =======  ======
 H       1.00e+00      +      +        +
He       8.51e-02      +      +        +
Li       1.26e-11
Be       2.51e-11
 B       3.55e-10
 C       3.31e-04      +      +        +
 N       8.32e-05      +      +
 O       6.76e-04      +      +        +
 F       3.63e-08
Ne       1.20e-04      +      +        +
Na       2.14e-06      +
Mg       3.80e-05      +      +
Al       2.95e-06      +
Si       3.55e-05      +      +
 P       2.82e-07
 S       2.14e-05      +      +
Cl       3.16e-07
Ar       2.51e-06      +
 K       1.32e-07
Ca       2.29e-06      +
Sc       1.48e-09
Ti       1.05e-07
 V       1.00e-08
Cr       4.68e-07
Mn       2.45e-07
Fe       3.16e-05      +      +
Co       8.32e-08
Ni       1.78e-06      +
Cu       1.62e-08
Zn       3.98e-08
=======  =========  ========  =======  ======

It shoud be noted that CHIANTI does not include a complete set of data for every ion of every element in this list.

