========
Tutorial
========

The ChiantiPy Approach
----------------------

Python is a modern, object-oriented programming language.  It provides a number of features that are useful to the programmer and the end user, such as, classes with methods (function-like) and attributes (data), and functions, among other things.  ChiantiPy has been constructed so that the primary means to calculate the spectral properties of ions and groups of ions is by way of Python classes.

More detailed information can be found in the **API** Reference.

ChiantiPy Classes
~~~~~~~~~~~~~~~~~

There are 7 basic classes that are provided by ChiantiPy.

  **ion**
  
    this class is very useful in itself and is the basic unit employed by all of the other classes
    
  **continuum**
    
    for calculating the free-free (*bremstrahlung*), free-bound (radiative recombination) continuum as well as the radiative loss rates due to these processes.
    
  **bunch**
  
    allows the user to specify a *bunch* of ions and to calculate the radiative properties of the selected ions in a group.  The ions can be specified by list of individual ions, list of elements, as well as by the minimum elemental abundance.  The properties of each ion are available, as with members of the **ion** class.  Among other things, the ratios of lines of different ions and elements can be calculated and then displayed.
    
  **spectrum**

    the **spectrum** class calculates the intensities of the line and continuum and them convolves the complete spectrum with a filter such as a Gaussian of specified width
    
    there are actually 3 *spectrum* classes.  Two of these all the use of multiple cpu cores to speed the calculation.  The basic **spectrum** class does not do multi-processes and is therefore compatible with most Python environments
    
  **mspectrum**
  
    **mspectrum** duplicates the calculations of the **spectrum** class but it employes the Python multiprocessing package in order to use multiple cpu cores to calculate the spectrum.  This class can be used in a basic Python shell, in a Python script, or in an IPython terminal.  It can not be used in either the Jupyter qtconsole or the Jupyter notebook.
    
  **ipymspectrum**
  
    this class employes the IPython ipyparallel module to provide access to multiple cpu cores.  It can only be used in the Jupyter qtconsole and the Jupyter notebook.
    
  **ioneq**
  
    this class allows the user to load and plot the ionization equilibria of a specfic element.  It can read the ionization equililbrium files in the CHIANTI $XUVTOP/ioneq directory.  Different ionization equilibria can be plotted against each other
    
    it is also possible to calculate the ionization equilibria of an individual element using the ionization and recombination rates in the CHIANTI database.  The results of this calculation can also be plotted agains existing calculations in the CHIANTI $XUVTOP/ioneq directory.
    
    
ChiantiPy Classes, Methods and Attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each of the ChiantiPy classes listed above has a number of methods for calculating various properties.  The results of these calculations are stored as attributes of the class that has been instantiatied (created).  In Python, all objects, which includes everything in Python, have introspection so that all methods and attributes can be discovered and used.  The IPython terminal and the Jupyter qtconsole both provide their own methods of easily displaying the methods and attributes.

Some methods in each class are more useful to the user than others.  For example, below the populate() method is demonstrated below.  However, it is generally not necessary for the user to use the populate() method.  Methods that need the ion population will make use that the Population attribute is available and, if not, use the populate() method to create it.

Below, the methods that are most likely of interest to users are listed below.  All of the available methods are presented and documented in the API section, a part of the ChiantiPy documentation.

  **ion**
  
    **popPlot()** method
    
        plots the level populations of the *top* (most highly populated) levels as a function of temperature and/or density.
  
    **gofnt()** method
    
        an interactive method that plots the most intense lines of an ion in a give wavelength range (*wvlRange*) and allows the user to select a line or several lines, that will be summed, and then plots the *GofT* function for the selected lines and saves these values in the **Gofnt** dictionary as an attribute of the ion object.
      
    **emiss()** method
    
      calculates the spectral line emissivities of the ion and saves these in the **Emiss** dictionary as an attribute of the ion object.
        
    **intensity()** method
    
      calculates the *intensity* of the specified ion as a function of temperature and density.  These properties are saved in the **Intensity** dictionary, available as an attribute of the ion.
      
    **intensityList()** method
    
      lists the spectral line intensities in a given wavelength range (*wvlRange*) in an interactive terminal or notebook
      
    **intensityPlot()** method
    
      plots the spectral line intensites in a given wavelength range (*wvlRange*) for the *top* most intensity lines.
      
    **intensityRatio()** method
    
      an interactive method that plots the most intense lines of an ion in a give wavelength range (*wvlRange*) and allows the user to select a pair of lines or a pair of lines to be summed and then plots the intensity ratio as a function of temperature and/or density.  The ratio is saved in the **IntensityRatio** dictionary as an attribute of the ion.
      
    **spectrum()** method
    
      calculates the spectrum of the ion as a function of wavelength.  The spectral line intensities are pass through a selectable filter to simulate the spectrometer line profile.  The spectrum is save in the **Spectrum** dictionary as an attribute of the ion.
    
    **ionizRate()** method
    
      calculates the ionization rate coefficient as a function of temperature.  The rate coefficient is save in the **IonizRate** dictionary as an attribute.  Uses the methods **diRate()** and **eaRate()** to first calculate the direct and excitation-ionization (ea) rate coefficients and sums them.
      
    **recombRate()** method
    
      calculates the recombination rate coefficient as a function of temperature.  The rate coefficient is save in the **RecombRate** dictionary as an attribute.    Uses the methods **rrRate()** and **drRate()** to first calculate the radiative recombination and dielectronic recombination rate coefficients and sums them.
    
  **ioneq**
  
    **load()** method
    
      reads a selected, existing ionization equilibrium calculation for a given element and saves it as a numpy array **Ioneq** as an attribute of the object.
      
    **calculate()** method
    
      calculates the ionization equilibrium of a selected element from the CHIANTI ionization and recombination rates for a specified temperature(s) and saves it as a numpy array **Ioneq** as an attribute of the object.
      
    **plot()** method
    
      plots the loaded or calculated ionization equilibrium.  Various parameters can be specified to plot only those aspects that are desired.  Can also plot an additional existing ionization equilibrium for comparison
      
  **bunch**
  
    the init method calculates the spectral line intensities for the selection of ions save the information in the **Intensity** dictionary as an attribute .  It does not calculate the continuum.
  
    beyond the init method, the **bunch** class inherits all of the following methods that are described under the **ion** class above
      
    **intensityList()**
      
    **intensityPlot()**
      
    **intensityRatio()**
      
    in addition, it inherits the following methods that are described under the **spectrum** class below
      
    **convolve()**
      
    **lineSpectrumPlot()**
      
    **spectrumPlot()**
      
  **spectrum**
  
    the init method calculates the spectral line intensities and the continuum due to the free-free (*bremstrahlung*), free-bound (radiative recombination), and two-photon processes.  The line intensities are convolved using the **convolve()** method (below).  The sum is saved in the **Spectrum** dictionary as an attribute.
  
    beyond the init method, the **spectrum** class also inherits the same methods as the **bunch** class including **intensityList()**, **intensityPlot**, and **intensityRatio**.
      
    **convolve()**
    
      convolves the line spectrum with specified filter from **ChiantiPy.tools.filters** using a specified width.
      
    **lineSpectrumPlot()**
    
      plots the convolved line spectrum as a function wavelength
      
    **spectrumPlot()**
    
      plots the spectrum calculated by the init method.  The summed (integrated) spectrum can be plotted or the spectrum for a specific temperature can be plotted.
      
  **mspectrum**
  
    the mspectrum behaves in the same way as the **spectrum** class except that it invokes the Python multiprocessing module so that the calculations are made using a specified number of cpu cores.  **mspectrum** can not be used in the Jupyter qtconsole or notebook.
    
  **ipymspectrum**
  
    the ipymspectrum behaves in the same way as the **spectrum** class except that it invokes the IPython ipyparallel module so that the calculations are made using a specified number of cpu cores.  **ipymspectrum** can only be used in the IPython terminal or the Jupyter qtconsole or notebook.


The ion class, basic properties
-------------------------------

Bring up a Python session, or better yet, an IPython session

::

  import ChiantiPy.core as ch
  fe14 = ch.ion('fe_14')

The fe14 object is instantiated with a number of methods and data.  Methods start with lowercase letters and attributes start with uppercase letters.  It is best not to simply import ion as there is a method with the same name in matplotlib.  A few examples:

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

If the fe14 **ion** object had be instantiated (created) with a temperature and an electron density, then many more attributes can be calculated.  For example, if the populate() method is used, it creates a dictionary attribute *Population*.  One thing to remember with Python is that capitalization matters.  


::

  import numpy as np
  t = 10.**(5.8+0.1*np.arange(11))
  dens = 1.e+9
  fe14 = ch.ion('fe_14')
  fe14.populate()
  fe14.Population.keys()
  >>['ci', 'protonDensity', 'popmat', 'eDensity', 'rec', 'population', 'temperature']

  fe14.Population['population'].shape
  >>(21, 739)

  '%10.2e'%(fe14.Temperature[10])
  >> '  2.00e+06'

  fe14.Population['population'][10,:5]
  >>array([ 8.71775703e-01, 1.27867444e-01, 4.91230626e-09, 4.29120495e-08, 1.35517895e-08])


gives the population of the first 5 of 739 levels of Fe XIV at a temperature of 2.00e+6

to be continued
