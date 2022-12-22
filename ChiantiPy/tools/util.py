"""
Utility functions

Notes
-----
Some of these functions can be replaced by roman numeral and periodic table lookup libraries.
some functions using os.walk can be replaced by os.path
"""
import os
import fnmatch
import numpy as np
from scipy import interpolate
from scipy.special import exp1
import ChiantiPy.tools.constants as const


def between(array,limits):
    """
    Find the indices of `array` corresponding to values in the range given by `limits`

    Parameters
    ----------

    array:  `list`, `ndarray`
        contains a list of elements such a wavelengths

    limits:  `list`, `ndarray`
        a 2 element array specifying the lower and upper limits of the array to be included


    """
    array = np.asarray(array)
    nlines = len(array)
#    hi = np.where(array >= limits[0],range(1,nlines+1),0)
#    lo = np.where(array <= limits[1],range(1,nlines+1),0)
    hi = np.where(array >= limits[0],list(range(1,nlines+1)),0)
    lo = np.where(array <= limits[1],list(range(1,nlines+1)),0)

    hilo = hi&lo
    out = [a -1  for a in hilo if a > 0]
    return out


def z2element(z):
    """
    Convert atomic number `z` to its elemental symbol.

    Parameters
    ----------
    z : `int`
        The atomic number/nuclear charge

    Returns
    -------

    element : `str`
            the abbreviated element name
    """
    if z-1 < len(const.El):
        element = const.El[z-1]
    else:
        element = 'unkn'
    return element


def spectroscopic2name(el,roman):
    """
    Convert from spectroscopic notation, e.g. Fe XI to 'fe_11'.

    Parameters
    ----------

    el : `str`
        Elemental symbol, e.g. Fe for iron
    roman : `str`
        Roman numeral spectroscopic symbol
    """
    elu = el.lower()
    romanu = roman.upper()
    idx = const.Ionstage.index(romanu)
    gname = elu+'_'+str(idx+1)
    return gname


def zion2name(z,ion, dielectronic=False):
    """
    Convert atomic number and ion number to generic name, e.g. (26,13) to 'fe_13'

    Parameters
    ----------

    z : `int`
        The atomic number/nuclear charge
    ion : `int`
        the ionization stage, 1 for the neutral, 2 for the first ionization stage, ...
    dielectronic : `bool`, optional

    Returns
    -------

    thisone : str
        the CHIANTI style ion name, such as 'fe_13'
    """
    if ion == 0:
        thisone = 0
    elif ion == z+2:
        thisone = 0
    elif (z-1 < len(const.El)) and (ion <= z+1):
        thisone = const.El[z-1]+'_'+str(ion)
        if dielectronic:
            thisone += 'd'
    else:
        # this should not actually happen
        thisone = 0
    return thisone


def zion2dir(z,ion, dielectronic=False, xuvtop=''):
    """
    Convert atomic number and ion number to CHIANTI database directory.

    Parameters
    ----------

    z : `int`
        The atomic number/nuclear charge
    ion : `int`
        the ionization stage, 1 for the neutral, 2 for the first ionization stage, ...
    dielectronic : `bool`, optional
    xuvtop : `str`, optional
        Set different CHIANTI database than the default

    Returns
    -------

    fname : str
        the CHIANTI directory where the file for the ion specified by z and ion are found
    """
    if xuvtop:
        dir = xuvtop
    else:
        dir = os.environ["XUVTOP"]
    if (z-1 < len(const.El)) and (ion <= z+1):
        thisel = const.El[z-1]
    else:
        thisel = ''
    if z-1 < len(const.El):
        thisone = const.El[z-1]+'_'+str(ion)
        if dielectronic:
            thisone += 'd'
    else:
        thisone = ''
    if thisel != '' :
        fname = os.path.join(dir,thisel,thisone)
    return fname


def zion2filename(z,ion, dielectronic=False, xuvtop=''):
    """
    Convert nuclear charge/atomic number and ion number to CHIANTI database filename.

    Parameters
    ----------
    z : `int`
        The atomic number/nuclear charge
    ion : `int`
        the ionization stage, 1 for neutrals, 2 for singly ionized
    dielectronic : `bool`, optional
        whether the ion is the simple dielectronic model
    xuvtop : `str`, optional
        the top directory of the CHIANTI database to be used

    Returns
    -------
    fname : str
        CHIANTI database filename
    """
    if xuvtop:
        dir = xuvtop
    else:
        dir = os.environ["XUVTOP"]
    if (z-1 < len(const.El)) and (ion <= z+1):
        thisel = const.El[z-1]
    else:
        thisel = ''
    if z-1 < len(const.El):
        thisone = const.El[z-1]+'_'+str(ion)
        if dielectronic:
            thisone += 'd'
    else:
        thisone = ''
    if thisel != '' :
        fname = os.path.join(dir,thisel,thisone,thisone)
    return fname


def zion2localFilename(z,ion, dielectronic=False):
    """
    Convert atomic number and ion number to generic file name with current directory at top.

    Parameters
    ----------

    z : `int`
        The atomic number/nuclear charge
    ion : `int`
        the ionization stage, 1 for neutrals, 2 for singly ionized
    dielectronic : `bool`, optional
    """
    dir = '.'
    if (z-1 < len(const.El)) and (ion <= z+1):
        thisel = const.El[z-1]
    else:
        thisel = ''
    if z-1 < len(const.El):
        thisone = const.El[z-1]+'_'+str(ion)
        if dielectronic:
            thisone += 'd'
    else:
        thisone = ''
    if thisel != '' :
        fname = os.path.join(dir,thisel,thisone,thisone)
    return fname

def zion2experimental(z, ion, dielectronic=False):
    """
    Convert atomic number and ion number to spectroscopic notation string

    Parameters
    ----------

    z : `int`
        The atomic number/nuclear charge
    ion : `int`
        the ionization stage, 0 for neutrals, 1 for singly ionized
    dielectronic : `bool`, optional

    Returns
    -------
    expt : `str`
        the experimental/laboratory notation for the ion, such as Fe 13+
    """
    if (z-1 < len(const.El)) and (ion <= z+1):
        expt = const.El[z-1].capitalize()+' ' + '%i'%(ion-1) + '+'
        if dielectronic:
            expt += ' d'
    else:  expt = ''
    return expt



def zion2spectroscopic(z,ion, dielectronic=False):
    """
    Convert atomic number and ion number to spectroscopic notation string

    Parameters
    ----------

    z : `int`
        The atomic number/nuclear charge
    ion : `int`
        the ionization stage, 1 for neutrals, 2 for singly ionized
    dielectronic : `bool`, optional

    Returns
    -------
    spect : `str`
        the spectroscopic notation for the ion, such as Fe XIV
    """
    if (z-1 < len(const.El)) and (ion <= z+1):
        spect = const.El[z-1].capitalize()+' '+const.Ionstage[ion-1]
        if dielectronic:
            spect += ' d'
    else:  spect = ''
    return spect


def convertName(name):
    """
    Convert ion name string (e.g. 'fe_13') to atomic number and ion number.

    Parameters
    ----------
    name : `str`
        the CHIANTI style name for an ion, such as fe_14
    Returns
    -------
    {'Z', 'Ion', 'Dielectronic', 'Element', 'higher', 'lower'} : `dict`
        `higher` and `lower` are the Chianti-style names for the higher and lower ionization stages, respectively.
    """
    s2 = name.split('_')
    els = s2[0].strip()
    Z = int(const.El.index(els)+1)
    ions = s2[1].strip()
    d = ions.find('d')
    if d > 0 :
        dielectronic = True
        ions = ions.replace('d','')
    else:
        dielectronic = False
    stage = int(ions)
    higher = zion2name(Z, stage+1)
    lower = zion2name(Z, stage-1)
    filename = zion2filename(Z, stage, dielectronic = dielectronic)
    iso = Z - stage + 1
    isoEl = const.El[iso - 1].capitalize()
    Zion = stage -1
    electrons = Z - Zion
    spectro = zion2spectroscopic(Z, stage)
    expt = zion2experimental(Z,  stage)
    return {'Z':Z,'Ion':stage, 'Zion':Zion, 'electrons':electrons, 'Dielectronic':dielectronic,
        'Element':els.capitalize(), 'higher':higher, 'lower':lower, 'filename':filename, 'iso':iso,
        'isoEl':isoEl, 'spectroscopic':spectro, 'experimental':expt}


def ion2filename(ions):
    """
    Convert ion name string to generic directory-file name.
    convertName has probably made this redundant

    Parameters
    ----------
    ions : str
        the CHIANTI style name for an ion, such as fe_14

    Returns
    -------
    fname:  str
        the full file name of the ion directory in the CHIANTI database
        assumes a top directory from the environmental variable XUVTOP

    Todo
    ----

    this duplicates what 'convertName' does
    """
    dir = os.environ["XUVTOP"]
    nameDict = convertName(ions)
    zion = nameDict['Z']
    el = z2element(zion)
    fname = os.path.join(dir,el,ions,ions)
    return fname


def el2z(els):
    """
    Convert elemental symbol to atomic number

    Parameters
    ----------

    els : str
        the abreviated element name

    Returns
    -------

    z : int
        the atomic number or nuclear charge of the element
    """
    z = const.El.index(els.lower())+1
    return z


def qrp(z,u):
    """
    Calculate :math:`Q_R^{\prime}(Z,u)`, where :math:`u=\epsilon/I` is the impact electron energy in threshold units, from Eq. 2.12 of [4]_.

    Parameters
    ----------
    z : `int`
        Atomic number
    u : array-like
        Impact electron energy in threshold units.

    Returns
    -------
    q : array-like
        1s ionization cross-section, :math:`Q_R^{\prime}(Z,u)`

    Notes
    -----
    Used for calculations of direct ionization cross-sections of the H and He sequences in `ChiantiPy.tools.io.twophotonHRead` and `ChiantiPy.tools.io.twophotonHeRead`, respectively.

    References
    ----------
    .. [4] Fontes, C. et al., 1999, PhRvA, `59, 1329 <http://adsabs.harvard.edu/abs/1999PhRvA..59.1329F>`_
    """
    #
    aa = 1.13  # aa stands for A in equ 2.12
    #
    if z >= 16 :
        # use Fontes Z=20, N=1 parameters
        dd = 3.70590
        c = -0.28394
        d = 1.95270
        cc = 0.20594
    else:
    # use Fontes Z=10, N=2 parameters
        dd = 3.82652
        c = -0.80414
        d = 2.32431
        cc = 0.14424
    #
    if z > 20:
        cc += ((z-20.)/50.5)**1.11
    #
    bu = u <= 1.
    q = np.ma.array(u, np.float64, mask=bu, fill_value=0.)
    #
    #
    q = (aa*np.ma.log(u) + dd*(1.-1./u)**2 + cc*u*(1.-1./u)**4 + (c/u+d/u**2)*(1.-1/u))/u
    #
    q.set_fill_value(0.)  # I don't know why this is necessary
    return q  #  .set_fill_value(0.)


def splomDescale(splom, energy):
    """
    Calculate the collision strength for excitation-autoionization as a function of energy.

    Parameters
    ----------
    energy : array-like
        In eV
    splom : `dict`
        Structure returned by `ChiantiPy.tools.io.splomRead`

    Returns
    -------
    omega : array-like
        Collision strength
    """
    #
    #
    nenergy = energy.size
    nsplom = len(splom['de'])
    # for these files, there are 5 spline points
    nspl = [splom['splom'].shape[0]]*splom['splom'].shape[1]
    if nenergy > 1:
        omega = np.zeros((nsplom,nenergy),np.float64)
    else:
        omega = np.zeros(nsplom,np.float64)
    #
    for isplom in range(nsplom):
        #
        dx = 1./(float(nspl[isplom])-1.)
        sxint = dx*np.arange(nspl[isplom])  # IDL sx
        sx1 = energy/(splom['de'][isplom]*const.ryd2Ev)  # IDL x_int
        good = sx1 >= 1.
        # make sure there are some valid energies above the threshold
        if good.sum():
            nbad = nenergy - good.sum()
            c_curr = splom['c'][isplom]
            #
            if splom['ttype'][isplom] == 1:
                sx = 1. - np.log(c_curr)/np.log(sx1[good] - 1. + c_curr)  # IDL sx_int
                y2 = interpolate.splrep(sxint,splom['splom'][:, isplom],s=0)  #allow smoothing,s=0)
                som = interpolate.splev(sx,y2,der=0)
                omega[isplom, nbad:] = som*np.log(sx1[good] -1. + np.exp(1.))
            #
            elif splom['ttype'][isplom] == 2:
                sx =(sx1[good] - 1.)/(sx1[good] -1. + c_curr)
                y2 = interpolate.splrep(sxint,splom['splom'][:, isplom],s=0)  #allow smoothing,s=0)
                som = interpolate.splev(sx,y2,der=0)
                omega[isplom, nbad:] = som
            #
            elif splom['ttype'][isplom] == 3:
                sx = (sx1[good] - 1.)/(sx1[good] -1. + c_curr)
                y2 = interpolate.splrep(sxint,splom['splom'][:, isplom],s=0)  #allow smoothing,s=0)
                som = interpolate.splev(sx,y2,der=0)
                omega[isplom, nbad:] = som/sx1[good]**2
            #
            elif splom['ttype'][isplom] == 4:
                sx = 1. - np.log(c_curr)/np.log(sx1[good] -1. + c_curr)
                y2 = interpolate.splrep(sxint,splom['splom'][:, isplom],s=0)  #allow smoothing,s=0)
                som = interpolate.splev(sx,y2,der=0)
                omega[isplom, nbad:] = som*np.log(sx1[good] -1. + c_curr)
            #
            #
            #
            elif splom['ttype'] > 4:
                print((' splom t_type ne 1,2,3,4 = %4i %4i %4i'%(splom['ttype'],splom['l1'],splom['l2'])))
        else:
            # there are no energies above the threshold
            pass
    #
    #
    omega = np.where(omega > 0.,omega,0.)
    #
    return omega


def dilute(radius):
    """
    Calculate the dilution factor.

    Parameters
    ----------
    radius : array-like
        Distance from the center of a star in units of the stellar radius.

    Notes
    -----
    If `radius` is less than 1.0, returns 0.
    """
    if radius >= 1.:
        d = 0.5*(1. - np.sqrt(1. - 1./radius**2))
    else:
        d = 0.
    return d


def listFiles(dir):
    """
    Walks the path and subdirectories to return a list of files.

    Parameters
    ----------

    dir : str
        the top directory to search
        subdirectories are also searched

    Returns
    -------

    listname:  list
        a list of files in dir and  subdirectories

    Notes
    -----
    This can be replaced by functions in `os.path`, as if 3.4, pathlib is probably better.
    It is not clear that this function is used anywhere in ChiantiPy
    """
    alist = os.walk(dir)
    listname = []
    for (dirpath,dirnames,filenames) in alist:
        if len(dirnames) == 0:
            for f in filenames:
                file = os.path.join(dirpath,f)
                if os.path.isfile(file):
                    listname.append(file)
        else:
            for f in filenames:
                file = os.path.join(dirpath,f)
                if os.path.isfile(file):
                    listname.append(file)
    return listname


def listRootNames(dir):
    """
    Walks the path and subdirectories to return a list of file root names.

    Notes
    -----
    This can be replaced by functions in `os.path`, as if 3.4, pathlib is probably better.
    Only seems to be used by
    """
    alist = os.walk(dir)
#    print(' getting file list')
    rootNames = []
    for (dirpath,dirnames,filenames) in alist:
        if len(dirnames) == 0:
            for f in filenames:
                file = os.path.join(dirpath,f)
                if os.path.isfile(file):
                    rootNames.append(os.path.splitext(f)[0])
        else:
            for f in filenames:
                file = os.path.join(dirpath,f)
                if os.path.isfile(file):
                    rootNames.append(os.path.splitext(f)[0])
    return rootNames

def findFileTypes(wpath, type ='*.txt',  verbose=False):
    """ to find all the files in wpath and below with file names matching fname
    """
    alist=sorted(os.walk(wpath))
    if verbose:
        print(' getting file list')
    listPath = []
    listFile = []
    fileName = []
    for (dirpath,dirnames,filenames) in alist:
        if len(dirnames) == 0:
            for f in filenames:
                if fnmatch.fnmatch(f, type):
                    file=os.path.join(dirpath,f)
                    if os.path.isfile(file):
                        listFile.append(file)
                        listPath.append(dirpath)
                        fileName.append(f)
        else:
            for f in filenames:
                if fnmatch.fnmatch(f, type):
                    file=os.path.join(dirpath,f)
                    if os.path.isfile(file):
                        listFile.append(file)
                        listPath.append(dirpath)
                        fileName.append(f)
    if verbose:
        nfiles=len(listFile)
        print((' nfiles = %i'%(nfiles)))
    return {'fullName':listFile, 'fileName':fileName}

def scale_bti(evin,crossin,f,ev1):
    """
    Apply ionization scaling of [7]_,[8]_, to energy and cross-section.

    Parameters
    ----------
    evin : float
        Energy - same units as ev1
    crossin : array-like
        Cross-section
    f : float -  the scale factor
    ev1 : float
        the ionization potential
        units - the same as evin

    Returns
    -------
    [bte,btx] : `list`
        Scaled energy and cross-section

    Notes
    -----
    This is the scaling used and discussed in the Dere (2007) calculation [1] of cross sections.  It was derived from similar scalings derived in reference [2]

    See Also
    --------
    descale_bti : Descale ionization energy and cross-section

    References
    ----------
    .. [7] Dere, K. P., 2007, A&A, `466, 771, <http://adsabs.harvard.edu/abs/2007A%26A...466..771D>`_
    .. [8] Burgess, A. and Tully, J. A., 1992, A&A, `254, 436 <http://adsabs.harvard.edu/abs/1992A%26A...254..436B>`_
    """
    u = evin/ev1
    bte = 1.-np.log(f)/np.log(u-1.+f)
    btx = u*crossin*(ev1**2)/(np.log(u)+1.)
    return [bte,btx]

def descale_bt(bte,btomega,f,ev1):
    """
    Apply excitation descaling of [3]_ to energy and collision strength

    Parameters
    ----------
    bte : array-like
        Scaled energy
    btomega : array-like
        Scaled collision strength
    f : array-like
    ev1 : array-like

    Returns
    -------
    [energy,omega] : `list`
        Descaled energy and collision strength

    Notes
    -----
    Need more details here. Not clear which equations are being used.

    See Also
    --------
    scale_bt : Apply scaling to excitation energy and cross-section

    References
    ----------
    .. [3] Burgess, A. and Tully, J. A., 1992, A&A, `254, 436 <http://adsabs.harvard.edu/abs/1992A%26A...254..436B>`_
    """
    u = 1.- f + np.exp(np.log(f)/(1.-bte))
    energy = u*ev1
    omega = (np.log(u) -1. + np.exp(1.))*btomega
    return [energy, omega]


def descale_bti(bte,btx,f,ev1):
    """
    Apply ionization descaling of [9]_ to energy and cross-sections of [7]_.

    Parameters
    ----------
    bte : array-like
        Scaled energy
    btx : array-like
        Scaled cross-section
    f : float
        Scaling parameter
    ev1 : float
        ionization potential - units determine the output energy units

    Returns
    -------
    [energy,cross] : `list`
        Descaled energy and cross-section

    Notes
    -----
    This is the scaling used and discussed in the Dere (2007) calculation [7]_ of cross sections.  It was derived from similar scalings provided by reference [2]

    See Also
    --------
    scale_bti : Descale ionization energy and cross-section

    References
    ----------
    .. [9] Burgess, A. and Tully, J. A., 1992, A&A, `254, 436 <http://adsabs.harvard.edu/abs/1992A%26A...254..436B>`_
    """
    u = 1.-f + np.exp(np.log(f)/(1. - bte))
    energy = u*ev1
    cross = (np.log(u)+1.)*btx/(u*ev1**2)
    return [energy,cross]


def scale_bt(evin,omega,f,ev1):
    """
    Apply excitation scaling of [8]_ to energy and collision strength.

    Parameters
    ----------
    evin : array-like
    omega : array-like
    f : array-like
    ev1 : array-like

    Returns
    -------
    [bte,btomega] : `list`
        Scaled energy and collision strength

    Notes
    -----
    Need more details here. Not clear which equations are being used.

    See Also
    --------
    descale_bt : Descale excitation energy and cross-section

    """
    u = evin/ev1
    bte = 1.-np.log(f)/np.log(u-1.+f)
    btomega = omega/(np.log(u)-1.+np.exp(1.))
    return [bte,btomega]

def descale_bti_rate(btTemperature, btRate, ip,  f=2.0):
    """
    Apply ionization descaling of [7]_, a Burgess-Tully type scaling to bt scaled ionization rates and
    temperatures.  The result is to return a temperature array and a ionization rate array.

    Parameters
    ----------
    btTemperature :  `array-like`
        the bt scaled temperature
    btRate : `array-like`
        the bt scaled ionization rate
    ip :  `float`
        the ionization potential in eV.
    f :  `float` (optional)
        the scaling parameter, 1.7 generally works well
    """
    temperature = (ip/const.boltzmannEv)*(np.exp(np.log(f/(1.-btTemperature))) - f)
    rate = (1./np.sqrt(btTemperature))*ip**(-1.5)*exp1(1./btTemperature)*btRate
    return {'temperature':temperature, 'rate':rate}

def scale_bt_rate(inDict, ip, f=2.0):
    """
    Apply ionization scaling of [7]_, a Burgess-Tully type scaling to ionization rates and
    temperatures. The result of the scaling is to return a scaled temperature between 0 and 1 and a
    slowly varying scaled rate as a function of scaled temperature. In addition, the scaled rates
    vary slowly along an iso-electronic sequence.

    Parameters
    ----------
    inDict : `dict`
        the input dictionary should have the following key pairs: `temperature`, array-like and
        `rate`, array-like
    ip :  `float`
        the ionization potential in eV.
    f :  `float` (optional)
        the scaling parameter, 1.7 generally works well

    Notes
    -----
    `btTemperature` and `btRate` keys are added to `inDict`
    """
    if ('temperature' and 'rate') in inDict.keys():
        rT = inDict['temperature']*const.boltzmannEv/ip
        btTemperature = 1. - np.log(f)/np.log(rT + f)
        btRate = np.sqrt(rT)*inDict['rate']*ip**1.5/(exp1(1./rT))
        inDict['btTemperature'] = btTemperature
        inDict['btRate'] = btRate
        inDict['ip'] = ip
    else:
        print(' input dict does not have the correct keys')
    return

def scale_classical(inDict, ip):
    """
    Apply the classical scaling to the input data

    Parameters
    ----------

    inDict: dictionary
        the input dictionary should have the following key pairs
            energy and cross
            or
            temperature and rate
    energy:  array-like
        energy values of the cross-section
    cross:  array-like
        a cross-section
    temperature:  array-like
    rate:  array-like
    ip:  float
        the ionization potential.  Typically in eV.

    Returns
        the following keys are added to inDict
    -------
    {'csEnergy', 'csCross', 'ip'} or {'csTemperature', 'csRate', 'ip'}
    """

    if ('energy' and 'cross') in inDict.keys():
        csEnergy = inDict['energy']/ip
        csCross = inDict['cross']*ip**2
        inDict['csEnergy'] = csEnergy
        inDict['csCross'] = csCross
        inDict['ip'] = ip
    elif ('temperature' and 'rate') in inDict.keys():
        csTemperature = inDict['temperature']/ip
        csRate = inDict['rate']*ip**2
        inDict['csTemperature'] = csTemperature
        inDict['csRate'] = csRate
        inDict['ip'] = ip
    else:
        print(' input dict does not have the correct keys')
    return

def units(defaults):
    """ to create a set of units compatible with ChiantiPy default values
    """
    # get xlabel for plotting
    if defaults['wavelength'] == 'angstrom':
        xlabel = 'Wavelength (\u212B)'
    elif defaults['wavelength'] == 'nm':
        xlabel = 'Wavelength (nm)'
    elif defaults['wavelength'] == 'ev':
        xlabel = 'Energy (eV)'
    elif defaults['wavelength'] == 'kev':
        xlabel = 'Energy (keV)'

    # get wavlength label for listing
    if defaults['wavelength'] == 'angstrom':
        listXlabel = 'Wvl (\u212B)'
    elif defaults['wavelength'] == 'nm':
        listXlabel = 'Wvl (nm)'
    elif defaults['wavelength'] == 'ev':
        listXlabel = 'Energy (eV)'
    elif defaults['wavelength'] == 'kev':
        listXlabel = 'Energy (keV)'
    #
    #  get Emiss Units

    if defaults['flux'] == 'energy':
        emissYlabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$ ion$^{-1}$'
    elif defaults['flux'] == 'photon':
        emissYlabel = r'photon cm$^{-2}$ s$^{-1}$ sr$^{-1}$ ion$^{-1}$'

    #  get Intensity Units

    if defaults['flux'] == 'energy':
        intYlabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$'
    elif defaults['flux'] == 'photon':
        intYlabel = r'photon cm$^{-2}$ s$^{-1}$ sr$^{-1}$'

    # get Spectrum Units

    if defaults['flux'] == 'energy':
        if defaults['wavelength'] == 'angstrom':
            spYlabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$' +  '\u212B' + '$^{-1}$'
        if defaults['wavelength'] == 'nm':
            spYlabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$ nm$^{-1}$'
        if defaults['wavelength'] == 'ev':
            spYlabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$ eV$^{-1}$'
        if defaults['wavelength'] == 'kev':
            spYlabel = r'erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$ keV$^{-1}$'
    elif defaults['flux'] == 'photon':
        if defaults['wavelength'] == 'angstrom':
            spYlabel = 'photon cm$^{-2}$ s$^{-1}$ sr$^{-1}$' +  '\u212B' + '$^{-1}$'
        if defaults['wavelength'] == 'nm':
            spYlabel = r'photon cm$^{-2}$ s$^{-1}$ sr$^{-1}$ nm$^{-1}$'
        if defaults['wavelength'] == 'ev':
            spYlabel = r'photon cm$^{-2}$ s$^{-1}$ sr$^{-1}$ eV$^{-1}$'
        if defaults['wavelength'] == 'kev':
            spYlabel = r'photon cm$^{-2}$ s$^{-1}$ sr$^{-1}$ keV$^{-1}$'

    # get radloss units

    radlossTlabel = 'Temperature (K)'
    radlossYlabel = 'erg  s$^{-1}$ cm$^{3}$'


    # get rate labels

    rateTlabel = 'Temperature (K)'
    rateYlabel = 's$^{-1}$ cm$^{3}$'

    labels = {}
    labels['xlabel'] = xlabel
    labels['listXlabel'] = listXlabel
    labels['emissYlabel'] = emissYlabel
    labels['intensityYlabel'] = intYlabel
    labels['spectrumYlabel'] = spYlabel
    labels['radlossTlabel'] = radlossTlabel
    labels['radlossYlabel'] = radlossYlabel
    labels['rateTlabel'] = rateTlabel
    labels['rateYlabel'] = rateYlabel

    return labels

