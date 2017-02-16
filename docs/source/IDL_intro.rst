===========================
Intro for CHIANTI IDL Users
===========================

ChiantiPy was not developed to be a clone of the CHIANTI IDL code.  The IDL code largely consists of functions that can be used to calculate or plot a variety of properties.  Structures are often used to carry the results of one function to be used by another.

ChiantiPy is object oriented.  For example, the building block of ChiantiPy is the **ion** class.  It carries with it all the methods that are needed as well as various calculated properties that are kept as attributes.  By following the Quick Start guide, you will become familiar with how ChiantiPy works.  Because one can always inquire, for example with `dir`, as to the methods and attributes of a an object, such as an ion, it is easy to remember what you might want to calculate next.  For example, you have created an object *myion*.  It is possible to then invoke

::

  myion.popPlot()
  
to plot the level populations of *myion*.  If you have not already calculated the levels populations, the **ion** class knows to calculate the level populations (*myion.populate()*) and save them for later use as the dictionary attribute *myion.Population* and then plot the specified level populations.  The level populations are then available as a 2 dimensional `numpy` array

::

  myion.Population['population']

Python and IPython provide many tools for examining any give object (everything in Python is an object of some sort)

::

  mg4 = ch.ion('mg_4',setup=0)
  
By setting the keyword *setup* to 0 or `False` the complete setup of the ion is not performed, but a certain amount of information is retrieved.  Since neither a temperature or electron density were specified, none of these attributes know anything related to these two quantitities.  The default value for setup is True so that the setup is performed with the specified temperatures and electron densities.  In general, this is the way you want to construct the mg_4 *ion*.

::

  for attr in vars(mg4):
      print(attr)
      
  IonStr
  FIP
  Z
  AbundanceName
  Ip
  Ion
  IoneqAll
  PDensity
  Abundance
  RadTemperature
  FileName
  RStar
  ProtonDensityRatio
  Dielectronic
  Defaults
  IoneqName
  Spectroscopic
  
The Python function `vars` retrieves the attributes of the mg4 object.

::

  print('%s'%(mg4.IoneqName))
  
  chianti
  
  print('the abundance file name is %s'%(mg4.AbundanceName))
  
  the abundance file name is sun_coronal_1992_feldman_ext
  
  print('the abundance of %s is %10.2e'%(mg4.Spectroscopic,mg4.Abundance))  
  
  the abundance of Mg IV is   1.41e-04
  
One can get a more complete description of the various attributes and methods of the mg4 object

::

  for one in dir(mg4):
      print(one)
      
  Abundance
  AbundanceName
  Defaults
  Dielectronic
  FIP
  FileName
  Ion
  IonStr
  IoneqAll
  IoneqName
  Ip
  PDensity
  ProtonDensityRatio
  RStar
  RadTemperature
  Spectroscopic
  Z
  __class__
  __delattr__
  __dict__
  __dir__
  __doc__
  __eq__
  __format__
  __ge__
  __getattribute__
  __gt__
  __hash__
  __init__
  __le__
  __lt__
  __module__
  __ne__
  __new__
  __reduce__
  __reduce_ex__
  __repr__
  __setattr__
  __sizeof__
  __str__
  __subclasshook__
  __weakref__
  boundBoundLoss
  cireclvlDescale
  convolve
  diCross
  diRate
  drRate
  drRateLvl
  eaCross
  eaDescale
  eaRate
  emiss
  emissList
  emissPlot
  emissRatio
  gofnt
  intensity
  intensityList
  intensityPlot
  intensityRatio
  intensityRatioInterpolate
  intensityRatioSave
  ionGate
  ioneqOne
  ionizCross
  ionizRate
  lineSpectrumPlot
  p2eRatio
  popPlot
  populate
  recombRate
  rrRate
  setup
  setupIonrec
  spectrum
  spectrumPlot
  twoPhoton
  twoPhotonEmiss
  twoPhotonLoss
  upsilonDescale
  upsilonDescaleSplups

First, at the top of the list are the attributes that were listed by the `vars` function.  Then come a number of methods starting with '__'.  These are generally not used and called *private* methods although nothing in Python is really private.  In IDL, just about everything is private.  After the *private* methods is a list of the methods provided by the mg4 *ion* class object.  These all start with a lower case letter to separate them from the attributes that start with an upper case letter (this is a ChiantiPy convention).

In the IPython and jupyter-qtconsole, typing 

::

  mg4.diCross(
  
and then hitting the tab key will bring up the doc-string for the diCross method, also found in the API reference in the documentation.  And then

::

  mg4.diCross()
  
calculates the direct ionization cross section of Mg IV for a set of energies above the ionization potential *Ip*.  The direct ionization cross sections are then provided in the mg4.DiCross dictionary.


Experience using the CHIANTI IDL package will provide the user with a background with what ChiantiPy can do.  However, the way to accomplish them are much easier but must be learned.  The best way to start is with the Quick Start guide and a book about Python.  Book suggestions are *Learning Python* by Mark Lutz and the handy *Python Pocket Reference*, also by Mark Lutz.  The first one is alway in reach and copies of the latter is on all of my computer desks.


