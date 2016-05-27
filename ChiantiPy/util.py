'''Utility functions, many for reading the CHIANTI database files.

Copyright 2009, 2010 Kenneth P. Dere

This software is distributed under the terms of the GNU General Public License
that is found in the LICENSE file


'''
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
import chianti
import chianti.constants as const
#import chianti.io as chio
#import chianti.gui as gui
#
#
def between(array,limits):
    '''
    returns an index array of elements of array where the values lie
    between the limits given as a 2 element list or tuple
    '''
    array=np.asarray(array)
    nlines=len(array)
#    hi=np.where(array >= limits[0],range(1,nlines+1),0)
#    lo=np.where(array <= limits[1],range(1,nlines+1),0)
    hi=np.where(array >= limits[0],list(range(1,nlines+1)),0)
    lo=np.where(array <= limits[1],list(range(1,nlines+1)),0)

    hilo=hi&lo
    out=[a -1  for a in hilo if a > 0]
    return out
    #
    # --------------------------------------------------
    #
def z2element(z):
    """ convert Z to element string """
    if z-1 < len(const.El):
        thisel=const.El[z-1]
    else:
        thisel=''
    return thisel
    #
    # -------------------------------------------------------------------------------------
    #
def spectroscopic2name(el,roman, dielectronic=False):
    """ 
    convert Z and ion to spectroscopic notation string 
    """
    elu = el.lower()
    romanu = roman.upper()
    idx =const.Ionstage.index(romanu)
    gname = elu+'_'+str(idx+1)
    return gname
    #
    # -------------------------------------------------------------------------------------
    #
def zion2name(z,ion, dielectronic=False):
    """
    convert Z, ion to generic name  26, 13 -> fe_13
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
    #
    # -------------------------------------------------------------------------------------
    #
def zion2dir(z,ion, dielectronic=False, xuvtop=0):
    """ convert Z to generic file name string """
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
    #
    # -------------------------------------------------------------------------------------
    #
def zion2filename(z,ion, dielectronic=False, xuvtop=0):
    """ convert Z to generic file name string """
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
    #
    # -------------------------------------------------------------------------------------
    #
def zion2localFilename(z,ion, dielectronic=False):
    """ convert Z to generic file name string with current directory at top"""
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
    """ convert ion name string to Z and Ion """
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
    convert ion string to generic file name string
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
    ''' qrp(Z,u)  u = E/IP
    calculate Qr-prime (equ. 2.12) of Fontes, Sampson and Zhang 1999'''
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
            elif ttype > 4:
                print((' splom t_type ne 1,2,3,4 = %4i %4i %4i'%(ttype,l1,l2)))
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
    #
    # ------------------------------------------------------------------------------
    #
def diCross(diParams, energy=0, verbose=0):
    '''
    Calculate the direct ionization cross section.
    diParams obtained by util.diRead with the following keys:
    ['info', 'ysplom', 'xsplom', 'btf', 'ev1', 'ref', 'eaev']
    Given as a function of the incident electron energy in eV
    returns a dictionary - {'energy':energy, 'cross':cross}
    '''
    iso=diParams['info']['iz'] - diParams['info']['ion'] + 1
    energy = np.array(energy, 'float64')
    if not energy.any():
        btenergy=0.1*np.arange(10)
        btenergy[0]=0.01
        dum=np.ones(len(btenergy))
        [energy, dum] = descale_bti(btenergy, dum, 2., diParams['ev1'][0])
        energy=np.asarray(energy, 'float64')
    #
    if iso == 1 and self.Z >= 6:
        #  hydrogenic sequence
        ryd=27.2113845/2.
        u=energy/self.Ip
        ev1ryd=self.Ip/ryd
        a0=0.5291772108e-8
        a_bohr=const.pi*a0**2   # area of bohr orbit
        if self.Z >= 20:
            ff = (140.+(self.Z/20.)**3.2)/141.
        else:
            ff = 1.
#        qr = util.qrp(self.Z,u)*ff
        qr = qrp(self.Z,u)*ff
        bb = 1.  # hydrogenic
        qh = bb*a_bohr*qr/ev1ryd**2
        diCross = {'energy':energy, 'cross':qh}
    elif iso == 2 and self.Z >= 10:
        #  use
        ryd=27.2113845/2.
        u=energy/self.Ip
        ev1ryd=self.Ip/ryd
        a0=0.5291772108e-8
        a_bohr=const.pi*a0**2   # area of bohr orbit
        if self.Z >= 20:
            ff=(140.+(self.Z/20.)**3.2)/141.
        else:
            ff=1.
#        qr=util.qrp(self.Z,u)*ff
        qr = qrp(self.Z,u)*ff
        bb = 2.  # helium-like
        qh=bb*a_bohr*qr/ev1ryd**2
        diCross={'energy':energy, 'cross':qh}
    else:
        cross=np.zeros(len(energy), 'Float64')

        for ifac in range(diParams['info']['nfac']):
            # prob. better to do this with masked arrays
            goode=energy > diParams['ev1'][ifac]
            if goode.sum() > 0:
                dum=np.ones(len(energy))
                btenergy, btdum = scale_bti(energy[goode],dum[goode], diParams['btf'][ifac], diParams['ev1'][ifac])
                # these interpolations were made with the scipy routine used here
                y2=interpolate.splrep(diParams['xsplom'][ifac], diParams['ysplom'][ifac], s=0)
                btcross=interpolate.splev(btenergy, y2, der=0)
                energy1, cross1 = descale_bti(btenergy, btcross, diParams['btf'][ifac], diParams['ev1'][ifac] )
                offset=len(energy)-goode.sum()
                if verbose:
                    pl.plot(diParams['xsplom'][ifac], diParams['ysplom'][ifac])
                    pl.plot(btenergy, btcross)
                if offset > 0:
                    seq=[np.zeros(offset, 'Float64'), cross1]
                    cross1=np.hstack(seq)
                cross+=cross1*1.e-14
        return {'energy':energy, 'cross':cross}
    #
    # ------------------------------------------------------------------------------
    #
def diCross1(diParams, energy=0, verbose=0):
    '''
    Calculate the direct ionization cross section.
    diParams obtained by util.diRead with the following keys:
    ['info', 'ysplom', 'xsplom', 'btf', 'ev1', 'ref', 'eaev']
    Given as a function of the incident electron energy in eV
    returns a dictionary - {'energy':energy, 'cross':cross}
    this version tests whether using the seq and hstack works
    so using a different approach
    '''
    iso=diParams['info']['iz'] - diParams['info']['ion'] + 1
    energy = np.array(energy, 'float64')
    if not energy.any():
        btenergy=0.1*np.arange(10)
        btenergy[0]=0.01
        dum=np.ones(len(btenergy))
        [energy, dum] = descale_bti(btenergy, dum, 2., diParams['ev1'][0])
        energy=np.asarray(energy, 'float64')
    #
    if iso == 1 and self.Z >= 6:
        #  hydrogenic sequence
        ryd=27.2113845/2.
        u=energy/self.Ip
        ev1ryd=self.Ip/ryd
        a0=0.5291772108e-8
        a_bohr=const.pi*a0**2   # area of bohr orbit
        if self.Z >= 20:
            ff = (140.+(self.Z/20.)**3.2)/141.
        else:
            ff = 1.
#        qr = util.qrp(self.Z,u)*ff
        qr = qrp(self.Z,u)*ff
        bb = 1.  # hydrogenic
        qh = bb*a_bohr*qr/ev1ryd**2
        diCross = {'energy':energy, 'cross':qh}
    elif iso == 2 and self.Z >= 10:
        #  use
        ryd=27.2113845/2.
        u=energy/self.Ip
        ev1ryd=self.Ip/ryd
        a0=0.5291772108e-8
        a_bohr=const.pi*a0**2   # area of bohr orbit
        if self.Z >= 20:
            ff=(140.+(self.Z/20.)**3.2)/141.
        else:
            ff=1.
#        qr=util.qrp(self.Z,u)*ff
        qr = qrp(self.Z,u)*ff
        bb=2.  # helium-like
        qh=bb*a_bohr*qr/ev1ryd**2
        diCross={'energy':energy, 'cross':qh}
    else:
        cross=np.zeros(len(energy), 'Float64')

        for ifac in range(diParams['info']['nfac']):
            # prob. better to do this with masked arrays
            goode=energy > diParams['ev1'][ifac]
            if goode.sum() > 0:
                dum=np.ones(len(energy))
                btenergy, btdum = scale_bti(energy[goode],dum[goode], diParams['btf'][ifac], diParams['ev1'][ifac])
                # these interpolations were made with the scipy routine used here
                y2=interpolate.splrep(diParams['xsplom'][ifac], diParams['ysplom'][ifac], s=0)
                btcross=interpolate.splev(btenergy, y2, der=0)
                energy1, cross1 = descale_bti(btenergy, btcross, diParams['btf'][ifac], diParams['ev1'][ifac] )
                offset=len(energy)-goode.sum()
                if verbose:
                    pl.plot(diParams['xsplom'][ifac], diParams['ysplom'][ifac])
                    pl.plot(btenergy, btcross)
#                if offset > 0:
#                    seq=[np.zeros(offset, 'Float64'), cross1]
#                    cross1=np.hstack(seq)
                cross[offset:]+=cross1*1.e-14
        return {'energy':energy, 'cross':cross}
    #
    # -------------------------------------------------------------------------------------
    #
def eaCross(diparams, easplom, elvlc, energy=None, verbose=False):
    '''
    Provide the excitation-autoionization cross section.

    Energy is given in eV.
    '''
    energy = np.asarray(energy, 'float64')
    if not energy.any():
        btenergy=0.1*np.arange(10)
        btenergy[0]=0.01
        dum=np.ones(len(btenergy))
        [energy, dum] = descale_bti(btenergy, dum, 2., min(easplom['deryd']))
        energy=np.asarray(energy, 'float64')
    #
    omega = splomDescale(easplom, energy)
    #
    #  need to replicate neaev

    if diparams['info']['neaev'] > 0:
        f1 = np.ones(omega.shape[0])
    else:
        f1 = diparams['info']['eaev']

    totalCross = np.zeros_like(energy)
    ntrans = omega.shape[0]
    for itrans in range(ntrans):
        lvl1 = easplom['lvl1'][itrans]
        mult = 2.*elvlc['j'][lvl1 - 1] + 1.
        cross = f1[itrans]*const.bohrCross*omega[itrans]/(mult.energy/const.ryd2Ev)
        totalCross += cross
    return {'energy':energy, 'cross':totalCross}
    #
    #-----------------------------------------------------------
    #
def eaDescale(easplups, temperature):
    """
    Calculates the effective collision strengths (upsilon)
    for excitation-autoionization as a function of temperature.
    a duplicate of ion.eaDescale()
    """
    #
    #  xt=kt/de
    #
    #  need to make sure elvl is >0, except for ground level
    #
    ntemp=temperature.size
    nsplups=len(easplups['de'])
    if ntemp > 1:
        ups=np.zeros((nsplups,ntemp),"Float64")
    else:
        ups=np.zeros(nsplups,"Float64")
    #
    for isplups in range(0,nsplups):
        l1=easplups["lvl1"][isplups]-1
        l2=easplups["lvl2"][isplups]-1
        ttype=easplups["ttype"][isplups]
        cups=easplups["cups"][isplups]
        nspl=easplups["nspl"][isplups]
        de=easplups["de"][isplups]
        dx=1./(float(nspl)-1.)
##                print easplups["easplups"][l1,l2]
        splups=easplups["splups"][isplups,0:nspl]
        kte=const.boltzmannEv*temperature/(const.ryd2Ev*de)
        #
        if ttype ==1:
            st=1.-np.log(cups)/np.log(kte+cups)
            xs=dx*np.arange(nspl)
            y2=interpolate.splrep(xs,splups,s=0)  #allow smoothing,s=0)
            sups=interpolate.splev(st,y2,der=0)
            ups[isplups]=sups*np.log(kte+np.exp(1.))
        #
        if ttype == 2:
            st=kte/(kte+cups)
            xs=dx*np.arange(nspl)
            y2=interpolate.splrep(xs,splups,s=0)
            sups=interpolate.splev(st,y2,der=0)
            ups[isplups]=sups
        #
        if ttype == 3:
            st=kte/(kte+cups)
            xs=dx*np.arange(nspl)
            y2=interpolate.splrep(xs,splups,s=0)
            sups=interpolate.splev(st,y2,der=0)
            ups[isplups]=sups/(kte+1.)
        #
        if ttype == 4:
            st=1.-np.log(cups)/np.log(kte+cups)
            xs=dx*np.arange(nspl)
            y2=interpolate.splrep(xs,splups,s=0)
            sups=interpolate.splev(st,y2,der=0)
            ups[isplups]=sups*np.log(kte+cups)
        #
        if ttype == 5:
            st=kte/(kte+cups)
            xs=dx*np.arange(nspl)
            y2=interpolate.splrep(xs,splups,s=0)  #allow smoothing,s=0)
            sups=interpolate.splev(st,y2,der=0)
            ups[isplups]=sups/(kte+0.)
        #
        #
        elif ttype > 5:  print((' t_type ne 1,2,3,4,5 = %5i %5i %5i'%(ttype,l1,l2)))
    #
    #
    ups=np.where(ups > 0.,ups,0.)
    #
#    easplups['ups'] = ups
    return ups
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
    returns [bte,btx]
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
