========
Tutorial
========


The ion class, basic properties
-------------------------------

Bring up a Python session, or better yet, an IPython session

::

  import chianti.core as ch
  fe14 = ch.ion('fe_14')

The fe14 object is instantiated with a number of methods and data.  Methods start with a lowercase letter and data start with uppercase.  It is best not to simply import ion as there is a method with the same name in matplotlib.  A few examples:

::

  fe14.IonStr
    >> 'fe_14'
  fe14.Spectroscopic
    >> 'Fe XIV'

CHIANTI and spectroscopic notation for the ion

::

  fe14.Z
    >> 26
  fe14.Ion
    >> 14
    
nuclear charge Z and the ionization stage (in spectroscopic notation) for the ion

::

  fe14.Ip
    >> 392.16196

this is the ionization potential in electron volts.

::

  fe14.FIP
    >> 7.9023801573028294

this is the first ionization potential (FIP) in electron volts - the ionization potential of the neutral (Fe I).

::

  fe14.Abundance
    >> 0.00012589265
  fe14.AbundanceName
    >> 'sun_photospheric_1998_grevesse'

this is the abundance of iron relative to hydrogen for the specified elemental abundance set.  For the ion class, the abundance can be specified by the *abuncance* keyword argument or the *abundanceName* keyword argument.  In the case the abundance is taken from default abundcance set.  The specified defaults can be examined by

::
	
  fe14.Defaults
    >>  {'abundfile': 'sun_photospheric_1998_grevesse', 'flux': 'energy', 'ioneqfile': 'chianti', 'wavelength': 'angstrom'}
 
the defaults can be specified by the user in the ~/.chiantirc/chiantrc file.  One is included in the distribution but it must be placed in ~/.chiantirc for it to be read.  If it is not found, a set of coded default values are used.

::

  fe14.Elvlc.keys()
  >>  ['ecmth', 'term', 'ref', 'pretty', 'spd', 'ecm', 'j', 'l', 'erydth', 'conf', 'lvl', 'spin', 'eryd', 'mult']

fe14.Elvlc is a dictionary that describes the energy levels of the Fe XIV ion.  The key 'ecm' provides the energies, relative to the ground level, in inverse cm.  The 'ref' key provides the references in the scientific literature where the data were provided.

::
	
  fe14.Elvlc['ref']
  >> ['%filename: fe_14.elvlc',
	%observed energy levels: Churilov S.S., Levashov V.E., 1993, Physica Scripta 48, 425,
	%observed energy levels: Redfors A., Litzen U., 1989, J.Opt.Soc.Am.B 6, #8, 1447,
	%theoretical energy levels: Storey P.J., Mason H.E., Young P.R., 2000, A&ASS 141, 28,
	%comment,
	Only level 16 does not have an observed energy. I have placed in,
	the third energy column a recommended value for the energy value of,
	this level, based on the theoretical and observed splittings of the,
	4F levels. It is this energy value which is used to compute the,
	wavelengths of transitions involving level 16 given in the .wgfa,
	file.,
	%produced as part of the Arcetri/Cambridge/GMU/NRL 'CHIANTI' atomic data base collaboration,
	%,
	%   P.R.Young Feb 99']

::

  fe14.Population.keys()
  >>['ci', 'protonDensity', 'popmat', 'eDensity', 'rec', 'population', 'temperature']

  fe14.Population['population'].shape
  >>(21, 739)

  '%10.2e'%(fe14.Temperature[10])
  >> '  2.00e+06'

  fe14.Population['population'][10,:5]
  >>array([  8.63535029e-01, 1.35791626e-01, 8.33350082e-09, 8.19101370e-08, 2.32292956e-08])


give the population of the first 5 of 739 levels of Fe XIV at a temperature of 2.00e+6

to be continued