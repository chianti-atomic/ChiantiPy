"""
Utility functions, many for reading the CHIANTI database files.

Notes
-----
Some of these functions can be replaced by roman numeral and periodic table lookup libraries.
"""

import os, fnmatch
#from types import *
#try:
#    # for Python 3 import
#    import configparser
#except ImportError:
#    # for Python 2 import
#    import ConfigParser as configparser
##from ConfigParser import *
import pickle
from datetime import date
import numpy as np
from scipy import interpolate
import ChiantiPy.tools.constants as const


def between(array,limits):
    """
    Find the indices of `array` corresponding to values in the range given by `limits`

    Parameters
    ----------
    array : `~numpy.ndarray`
    limits : `list` or `tuple` of length 2
    """
    array=np.asarray(array)
    nlines=len(array)
#    hi=np.where(array >= limits[0],range(1,nlines+1),0)
#    lo=np.where(array <= limits[1],range(1,nlines+1),0)
    hi=np.where(array >= limits[0],list(range(1,nlines+1)),0)
    lo=np.where(array <= limits[1],list(range(1,nlines+1)),0)

    hilo=hi&lo
    out=[a -1  for a in hilo if a > 0]
    return out


def z2element(z):
    """Convert atomic number `z` to its elemental symbol."""
    if z-1 < len(const.El):
        thisel=const.El[z-1]
    else:
        thisel=''
    return thisel


def spectroscopic2name(el,roman):
    """
    Convert from spectroscopic notation, e.g. Fe XI to 'fe_11'

    Parameters
    ----------
    el : `str`
        Elemental symbol, e.g. Fe for iron
    roman : `str`
        Roman numeral spectroscopic symbol
    """
    elu = el.lower()
    romanu = roman.upper()
    idx =const.Ionstage.index(romanu)
    gname = elu+'_'+str(idx+1)
    return gname


def zion2name(z,ion, dielectronic=False):
    """
    Convert atomic number and ion number to generic name, e.g. (26,13) to 'fe_13'

    Parameters
    ----------
    z : `int`
    ion : `int`
    dielectronic : `bool`, optional
    """
    if ion == 0:
        thisone = 0
    elif ion == z+2:
        thisone = 0
    elif (z-1 < len(const.El)) and (ion <= z+1):
        thisone=const.El[z-1]+'_'+str(ion)
        if dielectronic:
            thisone+='d'
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
    ion : `int`
    dielectronic : `bool`, optional
    xuvtop : `str`, optional
        Set different CHIANTI database than the default
    """
    if xuvtop:
        dir = xuvtop
    else:
        dir=os.environ["XUVTOP"]
    if (z-1 < len(const.El)) and (ion <= z+1):
        thisel=const.El[z-1]
    else:
        thisel=''
    if z-1 < len(const.El):
        thisone=const.El[z-1]+'_'+str(ion)
        if dielectronic:
            thisone+='d'
    else:
        thisone=''
    if thisel != '' :
        fname=os.path.join(dir,thisel,thisone)
    return fname


def zion2filename(z,ion, dielectronic=False, xuvtop=''):
    """
    Convert atomic number and ion number to CHIANTI database filename.

    Parameters
    ----------
    z : `int`
    ion : `int`
    dielectronic : `bool`, optional
    xuvtop : `str`, optional
        Set different CHIANTI database than the default
    """
    if xuvtop:
        dir = xuvtop
    else:
        dir=os.environ["XUVTOP"]
    if (z-1 < len(const.El)) and (ion <= z+1):
        thisel=const.El[z-1]
    else:
        thisel=''
    if z-1 < len(const.El):
        thisone=const.El[z-1]+'_'+str(ion)
        if dielectronic:
            thisone+='d'
    else:
        thisone=''
    if thisel != '' :
        fname=os.path.join(dir,thisel,thisone,thisone)
    return fname


def zion2localFilename(z,ion, dielectronic=False):
    """
    convert Z, ion to generic file name string with current directory at top
    """
    dir='.'
    if (z-1 < len(const.El)) and (ion <= z+1):
        thisel=const.El[z-1]
    else:
        thisel=''
    if z-1 < len(const.El):
        thisone=const.El[z-1]+'_'+str(ion)
        if dielectronic:
            thisone+='d'
    else:
        thisone=''
    if thisel != '' :
        fname=os.path.join(dir,thisel,thisone,thisone)
    return fname
    #
    # -------------------------------------------------------------------------------------
    #
def zion2spectroscopic(z,ion, dielectronic=False):
    """
    convert Z and ion to spectroscopic notation string
    """
    if (z-1 < len(const.El)) and (ion <= z+1):
        spect=const.El[z-1].capitalize()+' '+const.Ionstage[ion-1]
        if dielectronic:
            spect+=' d'
    else:  spect = ''
    return spect
    #
    # -------------------------------------------------------------------------------------
    #
def convertName(name):
    """
    convert ion name string to Z and Ion
    output is a dictionary with keys:
    'Z', 'Ion', 'Dielectronic', 'Element', 'higher', 'lower'
    'higher' is the Chianti-style name for the higher ionization stage
    'lower' is the Chianti-style name for the lower ionization stage
    """
    s2=name.split('_')
    els=s2[0].strip()
    i1=const.El.index(els)+1
    ions=s2[1].strip()
    d=ions.find('d')
    if d >0 :
        dielectronic=True
        ions=ions.replace('d','')
    else: dielectronic=False
    higher = zion2name(int(i1), int(ions)+1)
    lower = zion2name(int(i1), int(ions)-1)
    return {'Z':int(i1),'Ion':int(ions),'Dielectronic':dielectronic, 'Element':els, 'higher':higher, 'lower':lower}
    #
    # -------------------------------------------------------------------------------------
    #
def ion2filename(ions):
    """
    convert ion string to generic directory-file name string
    """
    dir=os.environ["XUVTOP"]
    zion=convertName(ions)
    el=z2element(zion['Z'])
    fname=os.path.join(dir,el,ions,ions)
    return fname
    #
    # -------------------------------------------------------------------------------------
    #
def el2z(els):
    """
    from an the name of the element (1-2 letter) return Z
    """
    z=const.El.index(els.lower())+1
    return z
    #
    # -------------------------------------------------------------------------------------
    #
def qrp(z,u):
    '''
    qrp(Z,u)  u = E/IP
    calculate Qr-prime (equ. 2.12) of Fontes, Sampson and Zhang 1999
    used for calculations of direct ionization cross sections of the H and He sequences
    '''
    #
    aa=1.13  # aa stands for A in equ 2.12
    #
    if z >= 16 :
        # use Fontes Z=20, N=1 parameters
        dd=3.70590
        c=-0.28394
        d=1.95270
        cc=0.20594
    else:
    # use Fontes Z=10, N=2 parameters
        dd=3.82652
        c=-0.80414
        d=2.32431
        cc=0.14424
    #
    if z > 20:
        cc+=((z-20.)/50.5)**1.11
    #
    bu=u <= 1.
    q=np.ma.array(u, 'Float64', mask=bu, fill_value=0.)
    #
    #
    q=(aa*np.ma.log(u) + dd*(1.-1./u)**2 + cc*u*(1.-1./u)**4 + (c/u+d/u**2)*(1.-1/u))/u
    #
    q.set_fill_value(0.)  # I don't know why this is necessary
    return q  #  .set_fill_value(0.)
    #
    #-----------------------------------------------------------
    #
def splomDescale(splom, energy):
    """
    Calculates the collision strength
    for excitation-autoionization as a function of energy.
    energy in eV
    splom is obtained by tools/io.splotRead
    """
    #
    #
    nenergy=energy.size
    nsplom=len(splom['deryd'])
    # for these files, there are 5 spline points
    nspl = 5
    if nenergy > 1:
        omega = np.zeros((nsplom,nenergy),"float64")
    else:
        omega = np.zeros(nsplom,"float64")
    #
    dx = 1./(float(nspl)-1.)
    sxint = dx*np.arange(nspl)
    for isplom in range(0,nsplom):
        #
        sx1 = energy/(splom['deryd'][isplom]*const.ryd2Ev)
        good = sx1 >= 1.
        # make sure there are some valid energies above the threshold
        if good.sum():
            nbad = nenergy - good.sum()
            c_curr = splom['c'][isplom]
            #
            if splom['ttype'][isplom] == 1:
                sx = 1. - np.log(c_curr)/np.log(sx1[good] - 1. + c_curr)
                y2 = interpolate.splrep(sxint,splom['splom'][:, isplom],s=0)  #allow smoothing,s=0)
                som = interpolate.splev(sx,y2,der=0)
                omega[isplom, nbad:] = som*np.log(sx -1. + np.exp(1.))
            #
            elif splom['ttype'][isplom] == 2:
                sx =(sx1[good] - 1.)/(sx1[good] -1. + c_curr)
                y2 = interpolate.splrep(sxint,splom['splom'][:, isplom],s=0)  #allow smoothing,s=0)
                som=interpolate.splev(sx,y2,der=0)
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
                som=interpolate.splev(sx,y2,der=0)
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
    omega=np.where(omega > 0.,omega,0.)
    #
    return omega
    #
    # ------------------------------------------------------------------------------
    #
def dilute(radius):
    '''
    to calculate the dilution factor as a function distance from
    the center of a star in units of the stellar radius
    a radius of less than 1.0 (incorrect) results in a dilution factor of 0.
    '''
    if radius >= 1.:
        d = 0.5*(1. - np.sqrt(1. - 1./radius**2))
    else:
        d = 0.
    return d
    #
    # -------------------------------------------------------------------------------------
    #
def listFiles(path):
    '''
    walks the path and subdirectories to return a list of files
    '''
    alist=os.walk(path)
#    print(' getting file list')
    listname=[]
    for (dirpath,dirnames,filenames) in alist:
        if len(dirnames) == 0:
            for f in filenames:
                file=os.path.join(dirpath,f)
                if os.path.isfile(file):
                    listname.append(file)
        else:
            for f in filenames:
                file=os.path.join(dirpath,f)
                if os.path.isfile(file):
                    listname.append(file)
    return listname
    #
    #-----------------------------------------------------------
    #
def scale_bti(evin,crossin,f,ev1):
    """
    apply BT ionization scaling to (energy, cross-section)
    returns BT scaled energy and cross-section [bte,btx]
    """
    u=evin/ev1
    bte=1.-np.log(f)/np.log(u-1.+f)
    btx=u*crossin*(ev1**2)/(np.log(u)+1.)
    return [bte,btx]
    #
    #-----------------------------------------------------------
    #
def descale_bti(bte,btx,f,ev1):
    """
    descale BT ionization scaling
    returns [energy,cross-section]
    """
    u=1.-f+np.exp(np.log(f)/(1.-bte))
    energy=u*ev1
    cross=(np.log(u)+1.)*btx/(u*ev1**2)
    return [energy,cross]
    #
    #-----------------------------------------------------------
    #
def descale_bt(bte,btomega,f,ev1):
    """
    descale BT excitation scaling
    returns [energy,collision strength]
    """
    u=1.-f+np.exp(np.log(f)/(1.-bte))
    energy=u*ev1
    omega=(np.log(u)-1.+np.exp(1.))*btomega
    return [energy,omega]
    #
    #-----------------------------------------------------------
    #
def scale_bt(evin,omega,f,ev1):
    """
    apply BT excitation scaling to (energy, collision strength)
    returns [bte,btomega]
    """
    u=evin/ev1
    bte=1.-np.log(f)/np.log(u-1.+f)
    btomega=omega/(np.log(u)-1.+np.exp(1.))
    return [bte,btomega]
