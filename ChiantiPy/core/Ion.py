"""
Ion class
"""
import os
import copy
import time

import numpy as np
from scipy.interpolate import splev, splrep
from scipy import special

import matplotlib.pyplot as plt

import ChiantiPy.tools.filters as chfilters
import ChiantiPy.tools.util as util
import ChiantiPy.tools.io as io
import ChiantiPy.tools.constants as const
import ChiantiPy.tools.data as chdata
import ChiantiPy.Gui as chGui
from ChiantiPy.base import ionTrails
from ChiantiPy.base import specTrails
from ChiantiPy.base import ioneqOne

#  upper level of two photon transitions for helium-like ions
heseqLvl2 = [-1,3,-1,-1,-1,3,6,6,-1,6,6,6,6,6,3,5,3,5,3,5,-1,5,-1,5,-1,5,-1,5,-1,5]


class ion(ioneqOne, ionTrails, specTrails):
    """
    The top level class for performing spectral calculations for an ion in the
    CHIANTI database.

    Parameters
    ----------
    ionStr : `str`
        CHIANTI notation for the given ion, e.g. 'fe_12' that corresponds to the `Fe XII` ion.
    temperature : `float` , `tuple`, `list`, `~numpy.ndarray`, optional
        Temperature array (Kelvin)
    eDensity : `float` , `tuple`, `list`, or `~numpy.ndarray`, optional
        Electron density array (:math:`\mathrm{cm^{-3}}` )
    pDensity : `float`, `tuple`, `list` or `~numpy.ndarray`, optional
        Proton density (:math:`\mathrm{cm}^{-3}` )
    radTemperature : `float` or `~numpy.ndarray`, optional
        Radiation black-body temperature (in Kelvin)
    rStar : `float` or `~numpy.ndarray`, optional
        Distance from the center of the star (in stellar radii)
    abundance : `float` or `str`, optional
        Elemental abundance relative to Hydrogen or name of CHIANTI abundance file
        to use, without the '.abund' suffix, e.g. 'sun_photospheric_1998_grevesse'.
    setup : `bool` or `str`, optional
        If True, run ion setup function. Otherwise, provide a limited number of
        attributes of the selected ion
    em : `float` or `~numpy.ndarray`, optional
        Emission Measure, for the line-of-sight emission measure
        (:math:`\mathrm{\int \, n_e \, n_H \, dl}`)
        (:math:`\mathrm{cm}^{-5}`.), for the volumetric emission measure
        :math:`\mathrm{\int \, n_e \, n_H \, dV}` (:math:`\mathrm{cm^{-3}}`).

    Attributes
    ----------
    IonStr : `str`
        Name of element plus ion, e.g. `fe_12` for Fe XII
    Z : `int`
        the nuclear charge, 26 for `fe_12`.
    Ion : `int`
        the ionization stage, 12 for `fe_12`.
    Dielectronic : `bool`
        true if the ion is a 'dielectronic' ion where the levels are populated
        by dielectronic recombination.
    Spectroscopic : `str`
        the spectroscopic notation for the ion, such as `Fe XII` for `fe_12`.
    Filename : `str`
        the complete name of the file `generic` filename in the CHIANTI
        database, such as `$XUVTOP/fe/fe_12/fe_12`.
    Ip : `float`
        the ionization potential of the ion
    FIP : `float`
        the first ionization potential of the element
    Defaults : `dict`
        these are specified by the software unless a *chiantirc* file is found
        in '$HOME/.chianti':

    Notes
    -----
    The keyword arguments temperature, eDensity, radTemperature, rStar, em must all be either a float or have the same dimension as the rest if specified as lists, tuples or arrays.

    The `Defaults` dict should have the following keys:

    * *abundfile*, the elemental abundance file, unless specified in
      *chiantirc* this defaults to *sun_photospheric_1998_grevesse*.
    * *ioneqfile*, the ionization equilibrium file name.  Unless specified
      in 'chiantirc' this is defaults to *chianti*.  Other choices are
      availble in $XUVTOP/ioneq
    * *wavelength*, the units of wavelength (Angstroms, nm, or kev), unless
      specified in the 'chiantirc' this is defaults to 'angstrom'.
    * *flux*, specified whether the line intensities are give in energy or
      photon fluxes, unless specified in the 'chiantirc' this is defaults to
      *energy*.
    * *gui*, specifies whether to use gui selection widgets (True) or to
      make selections on the command line (False).  Unless specified in the
      'chiantirc' this is defaults to *False*.

    """


    def __init__(self, ionStr, temperature=None, eDensity=None,
                pDensity='default', radTemperature=None, rStar=None,
                abundance=None, setup=True, em=None, verbose=0):

        self.IonStr = ionStr
        _tmp_convert_name = util.convertName(ionStr)
        self.Z = _tmp_convert_name['Z']
        self.Ion = _tmp_convert_name['Ion']
        self.Iso = _tmp_convert_name['iso']
        self.Dielectronic = _tmp_convert_name['Dielectronic']
        self.HigherName = _tmp_convert_name['higher']
        self.Spectroscopic = _tmp_convert_name['spectroscopic']
        self.FileName = _tmp_convert_name['filename']
        self.Defaults = chdata.Defaults
        self.Labels = util.units(chdata.Defaults)

        if abundance is not None:
            if isinstance(abundance,  float):
                self.Abundance = abundance
                # need the following for p2eRatio but will use default abundance
                self.AbundAll = chdata.AbundanceDefault['abundance']
            elif isinstance(abundance, str):
                ab = io.abundanceRead(abundance,  verbose=1)
                self.Abundance = ab['abundance'][self.Z-1]
                self.AbundanceName = abundance
                self.AbundAll = ab['abundance']

        else:
#            self.AbundanceName = self.Defaults['abundfile']
        #
#            ab = chdata.Abundance[self.AbundanceName]['abundance']
            self.Abundance = chdata.AbundanceDefault['abundance'][self.Z-1]
            self.AbundAll = chdata.AbundanceDefault['abundance']



        self.IoneqName = self.Defaults['ioneqfile']
        self.RadTemperature = radTemperature
        self.RStar = rStar

        #  ip in eV, but don't read for bare ions
        if self.Ion <= self.Z:
            self.Ip = chdata.Ip[self.Z-1, self.Ion-1-self.Dielectronic]
            self.FIP = chdata.Ip[self.Z-1, 0]
            if self.Dielectronic:
                self.UpperIp = chdata.Ip[self.Z-1, self.Ion-1]



        if setup and self.IonStr in chdata.MasterList:
            self.IoneqAll = chdata.IoneqAll
            #  this needs to go after setting temperature and reading ionization
            #  equilibria
            #  needs to know self.NTempDens first
            self.argCheck(temperature, eDensity, pDensity, em)
            self.ioneqOne()

            self.setup()
        else:
            if verbose and (self.IonStr not in chdata.MasterList):
                print(' ion %s not in masterlist, just various attributes, ionization, recombination rates'%(self.IonStr))
            if np.any(temperature) is not None:
                self.Temperature = np.atleast_1d(temperature)
                self.Ntemp = self.Temperature.size
                self.IoneqAll = chdata.IoneqAll
                self.argCheck(temperature, eDensity, pDensity, em)
                self.ioneqOne()
                self.setupIonrec()


        self.GrndLevels = chdata.GrndLevels[self.Iso]

    def setup(self, alternate_dir=None, verbose=False):
        """
        Setup various CHIANTI files for the ion including .wgfa, .elvlc, .scups,
        .psplups, .reclvl, .cilvl, and others.

        Parameters
        -----------

        alternate_dir : `str`
            directory cotaining the necessary files for a ChiantiPy ion; use to
            setup an ion with files not in the current CHIANTI directory
        verbose : `bool`

        Notes
        -----

        If ion is initiated with `setup=False`, call this method to do the
        setup at a later point.
        """
        if alternate_dir is not None:
            fileName = os.path.join(alternate_dir, self.IonStr)
            elvlcFileName = fileName+'.elvlc'
            wgfaFileName = fileName+'.wgfa'
        else:
            fileName = util.ion2filename(self.IonStr)
            elvlcFileName = None
            wgfaFileName = None
        if alternate_dir:
            if verbose:
                print(' using elvlcfile in alternate dir:  %s'%(elvlcFileName))
            self.Elvlc = io.elvlcRead('',filename=elvlcFileName)
            self.Wgfa = io.wgfaRead('',filename=wgfaFileName, elvlcname = elvlcFileName, total=True)
        else:
            self.Elvlc = io.elvlcRead(self.IonStr)
            self.Wgfa = io.wgfaRead(self.IonStr, total=True)
        self.Nlvls = len(self.Elvlc['lvl'])
        self.Nwgfa = len(self.Wgfa['lvl1'])
        nlvlWgfa = max(self.Wgfa['lvl2'])
        nlvlList = [nlvlWgfa]
        scupsfile = fileName + '.scups'
#        cilvlfile = fileName + '.cilvl'
#        reclvlfile = fileName + '.reclvl'
        rrlvlfile = fileName + '.rrlvl'
        autofile = fileName + '.auto'
        drParamsFile = fileName + '.drparams'
        rrParamsFile = fileName + '.rrparams'
        # read the scups/splups file
        if os.path.isfile(scupsfile):
            # happens the case of fe_3 and prob. a few others
            self.Scups = io.scupsRead(self.IonStr)
            self.Nscups = len(self.Scups['lvl1'])
            nlvlScups = max(self.Scups['lvl2'])
            nlvlList.append(nlvlScups)
        else:
            self.Nscups = 0
            nlvlScups = 0
        # ignoring these files for now, will just use rrlvl and auto for
        # calculation level resolved recombination rates
        # read cilvl file
        self.Ncilvl = 0
        #  .reclvl file may not exist
        self.Nreclvl = 0
#        if os.path.isfile(reclvlfile):
#            self.Reclvl = io.cireclvlRead(self.IonStr, filetype='reclvl')
#            self.Nreclvl = len(self.Reclvl['lvl1'])
#            nReclvl = max(self.Reclvl['lvl2'])
#            nlvlList.append(nReclvl)
#        else:
#            self.Nreclvl = 0
        if os.path.isfile(rrlvlfile):
            self.Rrlvl = io.cireclvlRead(self.IonStr, filetype='rrlvl')
            self.Nrrlvl = len(self.Rrlvl['lvl1'])
            nRrlvl = max(self.Rrlvl['lvl2'])
            self.Nrrlvl = nRrlvl
            nlvlList.append(nRrlvl)
        else:
            self.Nrrlvl = 0
        #  psplups file may not exist
        psplupsfile = fileName +'.psplups'
        if os.path.isfile(psplupsfile):
            self.Psplups = io.splupsRead(self.IonStr, filetype='psplups')
            self.Npsplups = len(self.Psplups["lvl1"])
        else:
            self.Npsplups = 0
        # drparams file may not exist
        if os.path.isfile(drParamsFile):
            self.DrParams = io.drRead(self.IonStr)
        if os.path.isfile(rrParamsFile):
            self.RrParams = io.rrRead(self.IonStr)

        #  .auto file may not exist
        if os.path.isfile(autofile):
            self.Auto = io.autoRead(self.IonStr, total=True)
            self.Nauto = len(self.Auto['lvl1'])
        else:
            self.Nauto = 0
        # need to determine the number of levels that can be populated
        nlvlElvlc = len(self.Elvlc['lvl'])
        #  elvlc file can have more levels than the rate level files
        self.Nlvls = min([nlvlElvlc, max(nlvlList)])

    def setupIonrec(self, alternate_dir=None, verbose=False):
        """
        Setup method for ion recombination and ionization rates.

        Notes
        ------
        Allows a bare-bones ion object to be setup up with just the ionization and recombination
        rates. For ions without a complete set of files - one that is not in the MasterList.
        """
        if hasattr(self, 'Temperature'):
            self.NTempDens = self.Temperature.size
        if alternate_dir:
            fileName = os.path.join(alternate_dir, self.IonStr)
        else:
            fileName = util.ion2filename(self.IonStr)
        elvlcname = fileName+'.elvlc'
        if os.path.isfile(elvlcname):
            self.Elvlc = io.elvlcRead('',elvlcname)
        else:
            zstuff = util.convertName(self.IonStr)
            if zstuff['Ion'] - zstuff['Z'] != 1:
                # don't expect one for the bare ion
                if verbose:
                    print(' Elvlc file missing for '+self.IonStr)
#            return
        cilvlFile = fileName +'.cilvl'
        if os.path.isfile(cilvlFile):
            self.Cilvl = io.cireclvlRead('',filename = fileName,
                                            filetype = 'cilvl')
            self.Ncilvl = len(self.Cilvl['lvl1'])
        else:
            self.Ncilvl = 0
        #  .reclvl file may not exist
        reclvlfile = fileName +'.reclvl'
        if os.path.isfile(reclvlfile):
            self.Reclvl = io.cireclvlRead('',filename=fileName, filetype='reclvl')
            self.Nreclvl = len(self.Reclvl['lvl1'])
        else:
            self.Nreclvl = 0

        drparamsFile = fileName +'.drparams'
        if os.path.isfile(drparamsFile):
           self.DrParams = io.drRead('',  filename=drparamsFile)

        rrparamsFile = fileName +'.rrparams'
        if os.path.isfile(rrparamsFile):
            self.RrParams = io.rrRead('', filename=rrparamsFile)

        diparamsFile = fileName +'.diparams'
        if os.path.isfile(diparamsFile):
            self.DiParams = io.diRead('',  filename=diparamsFile)

        easplomFile = fileName +'.easplom'
        if os.path.isfile(easplomFile):
            self.Easplom = io.splomRead('', ea=True, filename = easplomFile)

        easplupsFile = fileName +'.easplups'
        if os.path.isfile(easplupsFile):
            self.EaParams = io.eaRead('',  filename=easplupsFile)
            ups = self.eaDescale()
            self.EaParams['ups'] = ups

        if not hasattr(self, 'IoneqOne'):
            if verbose:
                print(' getting ioneqone')
            self.IoneqAll = chdata.IoneqAll
            self.ioneqOne()

    def diCross(self, energy=None, verbose=False):
        """
        Calculate the direct ionization cross section (cm$^2) as a function of the
        incident electron energy in eV, puts values into attribute DiCross

        Parameters
        ----------

        energy:  array-like
            incident electron energy in eV

        verbose:  bool, int
            with verbose set to True, printing is enabled

        Attributes
        ----------

            DiCross:  dict

                keys:  energy, cross
        """
        iso = self.Z - self.Ion + 1
        if energy is None:
            energy = self.Ip*10.**(0.025*np.arange(101))
        else:
            energy = np.asarray(energy, np.float64)
        #
        if iso == 1 and self.Z >= 6:
            #  hydrogenic sequence
            ryd = const.ryd2Ev
            u = energy/self.Ip
            ev1ryd = self.Ip/ryd
            a0 = 0.5291772108e-8
            a_bohr = const.pi*a0**2   # area of bohr orbit
            if self.Z >= 20:
                ff = (140.+(self.Z/20.)**3.2)/141.
            else:
                ff = 1.
            qr = util.qrp(self.Z,u)*ff
            bb = 1.  # hydrogenic
            qh = bb*a_bohr*qr/ev1ryd**2
            self.DiParams = {}
            self.DiParams['info'] = {'neaev':0}
            self.DiCross = {'energy':energy, 'cross':qh}
        elif iso == 2 and self.Z >= 10:
            #  use
            ryd = const.ryd2Ev
            u = energy/self.Ip
            ev1ryd = self.Ip/ryd
            a0 = 0.5291772108e-8
            a_bohr = const.pi*a0**2   # area of bohr orbit
            if self.Z >= 20:
                ff = (140.+(self.Z/20.)**3.2)/141.
            else:
                ff = 1.
            qr = util.qrp(self.Z,u)*ff
            bb = 2.  # helium-like
            qh = bb*a_bohr*qr/ev1ryd**2
            self.DiParams = {}
            self.DiParams['info'] = {'neaev':0}
            self.DiCross = {'energy':energy, 'cross':qh}
        else:
            if not hasattr(self, 'DiParams'):
                self.DiParams = io.diRead(self.IonStr)
            cross = np.zeros(len(energy), np.float64)
            for ifac in range(self.DiParams['info']['nfac']):
                # prob. better to do this with masked arrays
                goode = energy > self.DiParams['ev1'][ifac]
                if goode.sum() > 0:
                    dum = np.ones(len(energy))
                    btenergy, btdum = util.scale_bti(energy[goode],dum[goode],
                        self.DiParams['btf'][ifac], self.DiParams['ev1'][ifac])
                    # these interpolations were made with the scipy routine
                    # used here
                    y2 = splrep(self.DiParams['xsplom'][ifac],
                        self.DiParams['ysplom'][ifac], s=0)
                    btcross = splev(btenergy, y2, der=0)
                    energy1, cross1 = util.descale_bti(btenergy,
                                        btcross,self.DiParams['btf'][ifac],
                                        self.DiParams['ev1'][ifac] )
                    offset = len(energy)-goode.sum()
                    if verbose:
                        plt.plot(self.DiParams['xsplom'][ifac],
                                    self.DiParams['ysplom'][ifac])
                        plt.plot(btenergy, btcross)
                    if offset > 0:
                        seq = [np.zeros(offset, np.float64), cross1]
                        cross1 = np.hstack(seq)
                    cross += cross1*1.e-14
            self.DiCross = {'energy':energy, 'cross':cross}

    def diRate(self, ngl=20):
        """
        Calculate the direct ionization rate coefficient as a function of temperature (K)
        """

        if hasattr(self, 'Temperature'):
            temperature = self.Temperature
        else:
            print(' temperature is not defined in diRate')
            return
        # Gauss-Laguerre constants

        zgl = special.roots_laguerre(ngl)
        xgl = zgl[0]
        wgl = zgl[1]

#        xgl = const.xgl
#        wgl = const.wgl

        alpha = 5.287e+13
        tev = const.boltzmannEv*temperature

        nTempDens = self.NTempDens

        rate = np.zeros(nTempDens, np.float64)
        tset = set(temperature)
        if len(tset) == 1:
            t1 = temperature[0]
            tev1 = const.boltzmannEv*t1
            x0 = self.Ip/tev1 # Ip in eV
            beta = np.sqrt(const.boltzmann*t1)
            egl = self.Ip+xgl*tev1
            self.diCross(energy=egl)
            crossgl = self.DiCross['cross']
            term1 = wgl*xgl*crossgl
            term2 = wgl*crossgl
            newcross = alpha*beta*np.exp(-x0)*(term1.sum()+x0*term2.sum())
            rate1 = newcross
            for itemp in range(nTempDens):
                rate[itemp] = rate1
        else:
            for itemp in range(nTempDens):
                x0 = self.Ip/tev[itemp]  # Ip in eV
                beta = np.sqrt(const.boltzmann*temperature[itemp])
                egl = self.Ip+xgl*tev[itemp]
                self.diCross(energy=egl)
                crossgl = self.DiCross['cross']
                term1 = wgl*xgl*crossgl
                term2 = wgl*crossgl
                newcross = alpha*beta*np.exp(-x0)*(term1.sum()+x0*term2.sum())
                rate[itemp] = newcross
        self.DiRate = {'temperature':temperature, 'rate':rate,  'ngl':ngl}

    def eaDescale(self):
        """
        Calculates the effective collision strengths (upsilon)
        for excitation-autoionization as a function of temperature.
        """
        #  xt=kt/de
        #  need to make sure elvl is >0, except for ground level
        if hasattr(self, 'EaParams'):
            eaparams = self.EaParams
        else:
            self.EaParams = io.eaRead(self.IonStr)
            eaparams = self.EaParams

        if hasattr(self, 'Temperature'):
            temperature = self.Temperature
        else:
            print(' temperature is not defined in eaDescale')
            return
        ntemp = temperature.size
        nsplups = len(eaparams['de'])
        if ntemp > 1:
            ups = np.zeros((nsplups,ntemp),np.float64)
        else:
            ups = np.zeros(nsplups,np.float64)

        for isplups in range(0,nsplups):
            l1 = self.EaParams["lvl1"][isplups]-1
            l2 = self.EaParams["lvl2"][isplups]-1
            ttype = self.EaParams["ttype"][isplups]
            cups = self.EaParams["cups"][isplups]
            nspl = self.EaParams["nspl"][isplups]
            de = self.EaParams["de"][isplups]
            dx = 1./(float(nspl)-1.)
            splups = self.EaParams["splups"][isplups][0:nspl]
            kte = const.boltzmannEv*temperature/(const.ryd2Ev*de)
            if ttype == 1:
                st = 1.-np.log(cups)/np.log(kte+cups)
                xs = dx*np.arange(nspl)
                y2 = splrep(xs,splups,s=0)  #allow smoothing,s=0)
                sups = splev(st,y2,der=0)
                ups[isplups] = sups*np.log(kte+np.exp(1.))
            if ttype == 2:
                st = kte/(kte+cups)
                xs = dx*np.arange(nspl)
                y2 = splrep(xs,splups,s=0)
                sups = splev(st,y2,der=0)
                ups[isplups] = sups
            if ttype == 3:
                st = kte/(kte+cups)
                xs = dx*np.arange(nspl)
                y2 = splrep(xs,splups,s = 0)
                sups = splev(st,y2,der=0)
                ups[isplups] = sups/(kte+1.)
            if ttype == 4:
                st = 1.-np.log(cups)/np.log(kte+cups)
                xs = dx*np.arange(nspl)
                y2 = splrep(xs,splups,s=0)
                sups = splev(st,y2,der=0)
                ups[isplups] = sups*np.log(kte+cups)
            if ttype == 5:
                st = kte/(kte+cups)
                xs = dx*np.arange(nspl)
                y2 = splrep(xs,splups,s=0)  #allow smoothing,s=0)
                sups = splev(st,y2,der=0)
                ups[isplups] = sups/(kte+0.)
            elif ttype > 5:  print(' t_type ne 1,2,3,4,5 = %5i %5i %5i'%(ttype,l1,l2))

        ups = np.where(ups > 0.,ups,0.)
        self.EaParams['ups'] = ups
        return ups

    def eaCross(self, energy=None, verbose=False):
        """
        Provide the excitation-autoionization cross section.

        Energy is given in eV.
        """

        if not hasattr(self, 'DiParams'):
            self.DiParams = io.diRead(self.IonStr)
        if 'errorMessage' in self.DiParams.keys():
            #  this is a H-like or He-like ion
            # the excitation-autionization for this ion does not exist
            if verbose:
                print(' there is no EA cross section for this ion')
            self.EaCross = {'errorMessage':'there is no EA cross-section'}
            if energy is not None:
                self.EaCross['energy'] = energy
                self.EaCross['cross'] = np.zeros_like(energy)
            return
        elif self.DiParams['info']['neaev'] == 0:
            # the excitation-autionization for this ion does not exist
            if verbose:
                print(' there is no EA cross section for this ion')
            self.EaCross = {'errorMessage':'there is no EA cross-section'}
            if energy is not None:
                self.EaCross['energy'] = energy
                self.EaCross['cross'] = np.zeros_like(energy)
            return

        else:
            if hasattr(self, 'Easplom'):
                easplom = self.Easplom
            else:
                self.Easplom = io.splomRead(self.IonStr, ea=1)
                easplom = self.Easplom
            if energy is None:
                energy = self.Easplom['deryd'][0]*const.ryd2Ev*1.01*10.**(0.025*np.arange(101))
            # multiplicity of ground level already included
            #  splomDescale takes care of when energy < threshold
            omega = util.splomDescale(easplom, energy)
            #  need to replicate neaev
            ntrans = len(easplom['de'])
            eaev = self.DiParams['eaev']
            if len(eaev) == 1:
                for itrans in range(ntrans):
                    eaev.append(eaev[0])

            totalCross = np.zeros_like(energy)
#            ntrans = omega.shape[0]
            partialCross = np.zeros((ntrans, energy.size), np.float64)
            for itrans in range(ntrans):
                #  the collision strengths have already by divided by the
                #  statistical weight of the ground level 2j+1
                cross = eaev[itrans]*const.bohrCross*omega[itrans]/(energy/const.ryd2Ev)
                totalCross += cross
                partialCross[itrans] = cross
            self.EaCross = {'energy':energy, 'cross':totalCross,
                            'partial':partialCross}

    def eaRate(self):
        """
        Calculate the excitation-autoionization rate coefficient.
        """
        # get neaev from diparams file
        if not hasattr(self, 'DiParams'):
            self.DiParams = io.diRead(self.IonStr)
        if self.DiParams['info']['neaev'] == 0:
            #FIXME: raise an error here?
            # basically, this means that this method should not be invoked
            return
        else:
            if hasattr(self, 'Temperature'):
                temperature = self.Temperature
            else:
                btT = 0.1*np.arange(10)
                btT[0] = 0.01
                dum = np.ones(10, np.float64)
                [temperature, dum] = util.descale_bt(btT, dum, self.EaParams['cups'][0], self.DiParams['de'][0])
                self.Temperature = temperature
            if hasattr(self, 'EaParams'):
                eaparams = self.EaParams
            else:
                self.eaParams = io.eaRead(self.IonStr)
                self.eaDescale()
                eaparams = self.EaParams
            #  need to replicate neaev
            nups = len(eaparams['de'])
            tev = const.boltzmannEv*temperature
            ntemp = temperature.size
            partial = np.zeros((nups, ntemp), np.float64)
            earate = np.zeros(temperature.size, np.float64)
            eaev = self.DiParams['eaev']
            if len(eaev) == 1:
                for iups in range(nups):
                    eaev.append(eaev[0])

            for iups in range(nups):
                x0 = const.ryd2Ev*eaparams['de'][iups]/tev
                #  upsilon has already been divided by the statistical weight
                # of the ground level 2j+1
                earate1 = eaev[iups]*const.collision*eaparams['ups'][iups]*np.exp(-x0)/(np.sqrt(temperature))
                partial[iups] = earate1
                earate += earate1
            self.EaRate = {'rate':earate, 'temperature':temperature, 'partial':partial}

    def ionizCross(self, energy=None):
        """
        Provides the total ionization cross section.

        Notes
        -----
        uses `diCross`  and `eaCross`.
        """
        if energy is None:
            energy = self.Ip*1.01*10.**(0.025*np.arange(101))
        else:
            energy = np.asarray(energy, np.float64)

        if self.Z < self.Ion:
            self.IonizCross = {'cross':np.zeros_like(energy), 'energy':energy}
            return

        self.diCross(energy)
        if self.DiParams['info']['neaev'] == 0:
            ionizCross = self.DiCross['cross']
            self.EaCross = {'errorMessage':'there is no EA cross-section', 'energy':energy, 'cross':np.zeros_like(energy)}
        else:
            self.eaCross(energy)
            ionizCross = self.DiCross['cross'] + self.EaCross['cross']
        self.IonizCross = {'cross':ionizCross, 'energy':energy}

    def ionizRate(self, ngl=20):
        """
        Provides the total ionization rate.

        Calls diRate and eaRate.
        """
        if self.Z < self.Ion:
            self.IonizRate = {'rate':np.zeros_like(self.Temperature), 'temperature':self.Temperature}
            return
        self.diRate(ngl=ngl)
        self.eaRate()
        if self.DiParams['info']['neaev'] == 0:
            ionizrate = self.DiRate['rate']
        else:
            ionizrate = self.DiRate['rate']+self.EaRate['rate']
        self.IonizRate = {'rate':ionizrate,
                        'temperature':self.DiRate['temperature'],  'ngl':ngl}

    def rrRate(self):
        """
        Provide the radiative recombination rate coefficient as a function of temperature (K).
        """
        if hasattr(self, 'Temperature'):
            temperature = self.Temperature
        else:
            print(' temperature is not defined in rrRate')
            return

        rrparamsfile = util.ion2filename(self.IonStr) + '.rrparams'
        if hasattr(self, 'RrParams'):
            rrparams = self.RrParams
        elif os.path.isfile(rrparamsfile):
            self.RrParams = io.rrRead(self.IonStr)
            rrparams = self.RrParams
        else:
            self.RrRate = {'temperature':temperature, 'rate':np.zeros_like(temperature)}
            return

        if rrparams['rrtype'] == 1:
            a = rrparams['params'][3]
            b = rrparams['params'][4]
            t0 = rrparams['params'][5]
            t1 = rrparams['params'][6]
            rate = a/(np.sqrt(temperature/t0)*(1.+np.sqrt(temperature/t0))**(1.-b)*(1.+np.sqrt(temperature/t1))**(1.+b))
            self.RrRate = {'temperature':temperature, 'rate':rate}
        elif rrparams['rrtype'] == 2:
            a = rrparams['params'][3]
            b = rrparams['params'][4]
            t0 = rrparams['params'][5]
            t1 = rrparams['params'][6]
            c = rrparams['params'][7]
            t2 = rrparams['params'][8]
            b += c*np.exp(-t2/temperature)
            rate = a/(np.sqrt(temperature/t0)*(1.+np.sqrt(temperature/t0))**(1.-b)*(1.+np.sqrt(temperature/t1))**(1.+b))
            self.RrRate={'temperature':temperature, 'rate':rate}
        elif rrparams['rrtype'] == 3:
            a = rrparams['params'][2]
            b = rrparams['params'][3]
            rate = a/(temperature/1.e+4)**b
            self.RrRate = {'temperature':temperature, 'rate':rate}
        else:
            self.RrRate = {'temperature':temperature, 'rate':np.zeros_like(temperature)}

    def rrlvlDescale(self,  verbose=1):
        """
        Interpolate and extrapolate rrlvl rates.
        Used in level population calculations.
        """
        if hasattr(self, 'Temperature'):
            temperature = self.Temperature
        else:
            print(' temperature is not defined in rrlvlDescale')
            return {'errorMessage':' temperature is not defined in rrlvlDescale'}
        lvlfile = util.ion2filename(self.IonStr)+'.rrlvl'
        if hasattr(self, 'Rrlvl'):
            lvl = self.Rrlvl
        elif os.path.isfile(lvlfile):
            self.Rrlvl = io.cireclvlRead(self.IonStr, '.'+ 'rrlvl')
            lvl = self.Rrlvl
        else:
            self.RrlvlRate = {'rate':np.zeros_like(temperature)}
            return

        #  the rates and temperatures in rrlvl are not necessarily all the same
        nlvl = len(lvl['lvl1'])

#        rate = np.zeros(( nlvl, self.Ntemp), np.float64)
        rate = []
        for itrans in range(nlvl):
#            lvl2 = lvl['lvl2'][itrans]
            nTemp = lvl['ntemp'][itrans]
            y2 = splrep(np.log(lvl['temperature'][itrans, :nTemp-1]), np.log(lvl['rate'][itrans, :nTemp-1]))
            rr = np.exp(splev(np.log(self.Temperature),y2))
            rate.append(rr)
        rate = np.asarray(rate, np.float64)


        self.RrlvlRate = {'rate':rate, 'lvl1':lvl['lvl1'], 'lvl2':lvl['lvl2'], 'temperature':temperature, 'type':'rrlvl'}

    def drRate(self):
        """
        Provide the dielectronic recombination rate coefficient as a function of temperature (K).
        """
        if hasattr(self, 'Temperature'):
            temperature = self.Temperature
        else:
            print(' temperature is not defined in drRate')
            return {'errorMessage':' temperature is not defined in drRate'}
        drparamsfile = util.ion2filename(self.IonStr) + '.drparams'

        if hasattr(self, 'DrParams'):
            drparams = self.DrParams
        elif os.path.isfile(drparamsfile):
            self.DrParams = io.drRead(self.IonStr)
            drparams = self.DrParams
        else:
            self.DrRate = {'rate':np.zeros_like(temperature), 'temperature':temperature}
            return

        if drparams['drtype'] == 1:
            # badnell type
            drenergy = drparams['eparams']
            drcoef = drparams['cparams']
            gcoef = drenergy > 0.
            ncoef = gcoef.sum()
            rate = np.zeros(temperature.size, np.float64)
            for icoef in range(ncoef):
                rate += drcoef[icoef]*np.exp(-drenergy[icoef]/temperature)
            rate = rate/temperature**1.5
            self.DrRate = {'temperature':temperature, 'rate':rate}
        elif drparams['drtype'] == 2:
            # shull type
            params = drparams['params']
            adi = params[0]
            bdi = params[1]
            t0 = params[2]
            t1 = params[3]
            rate = adi*np.exp(-t0/temperature)*(1.+bdi*np.exp(-t1/temperature))/temperature**1.5
            self.DrRate = {'temperature':temperature, 'rate':rate}


    def drRateLvl(self, verbose=0):
        """
        to calculate the level resolved dielectronic rate from the higher ionization stage to
        the ion of interest
        rates are determined from autoionizing A-values
        the dictionary self.DrRateLvl contains
        rate = the dielectronic rate into an autoionizing level
        effRate = the dielectronic rate into an autoionizing level mutilplied by the branching
        ratio for a stabilizing transition
        totalRate = the sum of all the effRates
        """
        if not hasattr(self, 'Higher'):
            nameStuff = util.convertName(self.IonStr)
            z = nameStuff['Z']
            stage = nameStuff['Ion']
            higherStr = util.zion2name(z, stage+1)
            self.Higher = ion(higherStr, self.Temperature, self.EDensity)

#        coef2 = (const.planck)**3/(2.*const.pi*const.emass*const.boltzmann*self.Temperature)**1.5
        coef3 = (const.bohrCross*4.*const.ryd2erg/(const.boltzmann*self.Temperature))**1.5
        coef = coef3
        nt = self.Temperature.size
        allRate = []
        effRate = []
        de = []
        erg = []
        ipErg = self.Ip*const.ev2Erg
        lvl = []
        branch = []
        dekt = []
        totalRate = np.zeros(nt, np.float64)
        for i, avalue in enumerate(self.Auto['avalue']):
            elvl1idx = self.Elvlc['lvl'].index(self.Auto['lvl1'][i])
            if elvl1idx == 1:
                elvl2idx = self.Elvlc['lvl'].index(self.Auto['lvl2'][i])
                gUpper = float(self.Elvlc['mult'][elvl2idx])
                gLower = float(self.Higher.Elvlc['mult'][elvl1idx])
                ecm2 = self.Elvlc['ecm'][elvl2idx]
                if ecm2 < 0.:
                    ecm2 = self.Elvlc['ecmth'][elvl2idx]
                de1 = ecm2*const.invCm2Erg - self.Ip*const.ev2Erg
                erg.append(ecm2*const.invCm2Erg)
                de.append(de1)
                dekt1 = de1/(const.boltzmann*self.Temperature)
                dekt.append(dekt1)
                expkt = np.exp(-dekt1)
                rate = coef*gUpper*expkt*avalue/(2.*gLower)
                branch1 = self.Wgfa['avalueLvl'][elvl2idx]/(avalue + self.Wgfa['avalueLvl'][elvl2idx])
                branch.append(branch1)
                lvl.append(self.Auto['lvl2'][i])
                allRate.append(rate)
                effRate.append(rate*branch1)
                totalRate += rate*branch1
        self.DrRateLvl = {'rate':allRate, 'effRate':effRate, 'totalRate':totalRate,  'de':de, 'avalue':self.Auto['avalue'], 'lvl':lvl, 'branch':branch, 'dekt':dekt, 'erg':erg, 'ipErg':ipErg}

    def recombRate(self):
        """
        Provides the total recombination rate coefficient.

        Calls `drRate` and `rrRate`
        """

        if hasattr(self, 'Temperature'):
            temperature = self.Temperature
        else:
            print(' temperature is not defined in recombRate ion = %s'%(self.IonStr))
            self.RecombRate = {'errorMessage':' temperature is not defined in recombRate'}
        if self.Ion == 1:
            self.RecombRate = {'rate':np.zeros_like(temperature), 'temperature':temperature}
            return
        self.rrRate()
        self.drRate()
        rate = np.zeros_like(temperature)
        if hasattr(self, 'DrRate'):
            rate += self.DrRate['rate']
        if hasattr(self, 'RrRate'):
            rate += self.RrRate['rate']
        self.RecombRate = {'rate':rate, 'temperature':temperature}

    def upsilonDescale(self, prot=0):
        """
        Provides the temperatures and effective collision strengths (upsilons)
        set prot for proton rates
        otherwise, ce will be set for electron collision rates
        uses the new format "scups" files
        """

        if prot:
            ce = 0
            if hasattr(self, 'Psplups'):
                nscups = len(self.Psplups["lvl1"])
            else:
                self.Psplups = io.splupsRead(self.IonStr,filetype='psplups')
                if self.Psplups is None:
                    self.PUpsilon = None
                    return
                else:
                    nscups = len(self.Cilvl["lvl1"])
        else:
            ce=1
            if hasattr(self, 'Scups'):
                nscups = len(self.Scups["lvl1"])
            else:
                self.Scups = io.scupsRead(self.IonStr)
                if not self.Scups['status']:
                    self.Upsilon = None
                    return
                else:
                    nscups = len(self.Scups["lvl1"])

        if hasattr(self, 'Temperature'):
            temperature = self.Temperature
        else:
            print(' Temperature undefined in upsilonDescale')
            return {'errorMessage':' Temperature undefined in upsilonDescale'}

        if not hasattr(self, 'Elvlc'):
            self.elvlcRead()

        eryd = np.asarray(self.Elvlc["eryd"])
        erydth = np.asarray(self.Elvlc["erydth"])
        elvlc = np.where(eryd >= 0.,eryd,erydth)
        temp = np.asarray(temperature)
        ntemp = temp.size
        if ntemp > 1:
            ups = np.zeros((nscups,ntemp),np.float64)
            exRate = np.zeros((nscups,ntemp),np.float64)
            dexRate = np.zeros((nscups,ntemp),np.float64)
        else:
            ups = np.zeros(nscups,np.float64)
            exRate = np.zeros((nscups,ntemp),np.float64)
            dexRate = np.zeros((nscups,ntemp),np.float64)
        deAll = []
        for iscups in range(nscups):
            if prot:
                # for proton rates
                l1 = self.Psplups["lvl1"][iscups]-1
                l1idx = self.Elvlc['lvl'].index(self.Psplups['lvl1'][iscups])
                l2 = self.Psplups["lvl2"][iscups]-1
                l2idx = self.Elvlc['lvl'].index(self.Psplups['lvl2'][iscups])
                ttype = self.Psplups["ttype"][iscups]
                cups = self.Psplups["cups"][iscups]
                nspl = self.Psplups["nspl"][iscups]
                dx = 1./(float(nspl)-1.)
                xs = dx*np.arange(nspl)
                scups = self.Psplups["splups"][iscups]
                de = elvlc[l2idx]-elvlc[l1idx]
                kte = const.boltzmann*temp/(de*const.ryd2erg)
            else:
                # electron collisional excitation or dielectronic excitation
                l1 = self.Scups["lvl1"][iscups]-1
                l1idx = self.Elvlc['lvl'].index(self.Scups['lvl1'][iscups])
                l2 = self.Scups["lvl2"][iscups]-1
                l2idx = self.Elvlc['lvl'].index(self.Scups['lvl2'][iscups])
                ttype = self.Scups["ttype"][iscups]
                cups = self.Scups["cups"][iscups]
                nspl = self.Scups["ntemp"][iscups]
                xs = self.Scups['btemp'][iscups]
                scups = self.Scups["bscups"][iscups]
                de = self.Scups['de'][iscups]
                kte = const.boltzmann*temp/(de*const.ryd2erg)
            der=0
            if ttype == 1:
                st = 1.-np.log(cups)/np.log(kte+cups)
                y2 = splrep(xs,scups,s=0)
                sups = splev(st,y2,der=der)
                ups[iscups] = sups*np.log(kte+np.exp(1.))
            if ttype == 2:
                st = kte/(kte+cups)
                y2 = splrep(xs,scups,s = 0)
                sups = splev(st,y2,der = der)
                ups[iscups] = sups
            if ttype == 3:
                st = kte/(kte+cups)
                y2 = splrep(xs,scups,s = 0)
                sups = splev(st,y2,der=der)
                ups[iscups] = sups/(kte+1.)
            if ttype == 4:
                st = 1.-np.log(cups)/np.log(kte+cups)
                y2 = splrep(xs,scups,s = 0)
                sups = splev(st,y2,der=der)
                ups[iscups] = sups*np.log(kte+cups)
            if ttype == 5:
                # dielectronic rates
                st = kte/(kte+cups)
                y2 = splrep(xs,scups,s=0)
                sups = splev(st,y2,der=der)
                ups[iscups] = sups/(kte+0.)
            #  descale proton values
            if ttype == 6:
                st = kte/(kte+cups)
                y2 = splrep(xs,scups,s=0)
                sups = splev(st,y2,der=der)
                ups[iscups] = 10.**sups
            elif ttype > 6:  print(' t_type ne 1,2,3,4,5 = %5i %5i %5i '%(ttype,l1,l2))

            if ce:
                if self.Dielectronic:
                    # the dielectronic ions will eventually be discontinued
                    de = np.abs((elvlc[l2idx] - self.UpperIp/const.ryd2Ev) - elvlc[l1idx])
                else:
                    de = np.abs(elvlc[l2idx] - elvlc[l1idx])
                deAll.append(de)
                ekt = (de*const.ryd2erg)/(const.boltzmann*temp)
                fmult1 = float(self.Elvlc["mult"][l1idx])
                fmult2 = float(self.Elvlc["mult"][l2idx])
                dexRate[iscups] = const.collision*ups[iscups]/(fmult2*np.sqrt(temp))
                exRate[iscups] = const.collision*ups[iscups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
            elif prot:
                de = np.abs(elvlc[l2idx]- elvlc[l1idx])
                ekt = (de*1.57888e+5)/temp
                fmult1 = float(self.Elvlc["mult"][l1idx])
                fmult2 = float(self.Elvlc["mult"][l2idx])
                dexRate[iscups] = const.collision*ups[iscups]/(fmult2*np.sqrt(temp))
                exRate[iscups] = const.collision*ups[iscups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))

        ups=np.where(ups > 0.,ups,0.)
        if prot == 1:
            self.PUpsilon = {'upsilon':ups, 'temperature':temperature,
                                'exRate':exRate, 'dexRate':dexRate}
        else:
            self.Upsilon = {'upsilon':ups, 'temperature':temperature,
                            'exRate':exRate, 'dexRate':dexRate, 'de':deAll}

    def spectrum(self, wavelength, filter=(chfilters.gaussianR,1000.), label=0, allLines=1):
        """
        Calculates the line emission spectrum for the specified ion.

        Convolves the results of intensity to make them look like an observed spectrum
        the default filter is the gaussianR filter with a resolving power of 1000.  Other choices
        include chianti.filters.box and chianti.filters.gaussian.  When using the box filter,
        the width should equal the wavelength interval to keep the units of the continuum and line
        spectrum the same.

        includes ionization equilibrium and elemental abundances

        can be called multiple times to use different filters and widths
        uses label to keep the separate applications of spectrum sorted by the label
        for example, do .spectrum( .... label = 'test1')
        and do          .spectrum( ....  label = 'test2')
        then will get self.Spectrum.keys() = test1, test2 and
        self.Spectrum['test1'] = {'intensity':aspectrum,  'wvl':wavelength, 'filter':useFilter.__name__,
            'filterWidth':useFactor}

        Notes
        ------
        scipy.ndimage.filters also includes a range of filters.
        """
        self.Wavelength = wavelength
        if not hasattr(self, 'WvlRange'):
            self.WvlRange = [wavelength[0],  wavelength[-1]]

        useFilter = filter[0]
        useFactor = filter[1]
        wvlRange = [wavelength.min(), wavelength.max()]
        if hasattr(self, 'Intensity'):
            intensity = self.Intensity
        else:
            self.intensity(allLines=allLines)
            intensity = self.Intensity
        #  if intensity had been called with em, then the intensities are
        # already multiply by em
        if 'errorMessage' in intensity.keys():
            errorMessage = intensity['errorMessage']
        else:
            errorMessage = None
        if hasattr(self, 'Em'):
            em = self.Em
            useEm = 0

        # unicode character for angstrom is \u212B
#        if self.Em.max() == 1.:
#            ylabel = 'erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$ \u212B$^{-1}$ ($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'
#        else:
#            ylabel = 'erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$ \u212B$^{-1}$'

        xlabel = self.Labels['xlabel']
        ylabel = self.Labels['spectrumYlabel']

        if np.array_equal(self.Em, np.ones_like(self.Em)):
            ylabel += '($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'

#        if self.Defaults['wavelength'] == 'angstrom':
#            xlabel = 'Wavelength (\u212B)'
#        else:
#            xlabel = 'Wavelength ('+self.Defaults['wavelength'] +')'


        aspectrum = np.zeros((self.NTempDens, wavelength.size), np.float64)
        if not 'errorMessage' in self.Intensity.keys():
            idx = util.between(self.Intensity['wvl'], wvlRange)
            if len(idx) == 0:
#                    print(' no lines in wavelength range %12.2f - %12.2f'%(wavelength.min(), wavelength.max()))
#                    self.Spectrum = {'errorMessage':' no lines in wavelength range %12.2f - %12.2f'%(wavelength.min(), wavelength.max())}
                errorMessage =  '%s no lines in wavelength range %12.2f - %12.2f'%(self.IonStr, wavelength.min(), wavelength.max())
            else:
                for itemp in range(self.NTempDens):
                    for iwvl in idx:
                        wvlCalc = self.Intensity['wvl'][iwvl]
                        aspectrum[itemp] += useFilter(wavelength, wvlCalc, factor=useFactor)*self.Intensity['intensity'][itemp, iwvl]
                    if useEm:
                        aspectrum[itemp] *= em[itemp]
        integrated = aspectrum.sum(axis=0)

        if type(label) == type(''):
            if hasattr(self, 'Spectrum'):
                self.Spectrum[label] = {'intensity':aspectrum.squeeze(), 'integrated':integrated,
                    'wavelength':wavelength, 'filter':useFilter, 'filterWidth':useFactor,
                    'allLines':allLines, 'em':em, 'xlabel':xlabel, 'ylabel':ylabel}
                if errorMessage is not None:
                    self.Spectrum[label]['errorMessage'] = errorMessage
            else:
                self.Spectrum = {label:{'intensity':aspectrum.squeeze(), 'integrated':integrated,
                    'wavelength':wavelength, 'filter':useFilter, 'filterWidth':useFactor,
                    'allLines':allLines, 'em':em, 'xlabel':xlabel, 'ylabel':ylabel}}
                if errorMessage is not None:
                    self.Spectrum[label]['errorMessage'] = errorMessage
        else:
            self.Spectrum = {'intensity':aspectrum.squeeze(), 'integrated':integrated,
                'wavelength':wavelength, 'filter':useFilter, 'filterWidth':useFactor,
                'allLines':allLines, 'em':em, 'xlabel':xlabel, 'ylabel':ylabel}
            if errorMessage is not None:
                self.Spectrum['errorMessage'] = errorMessage



    def populate(self, popCorrect=True, verbose=False):
        """
        Calculate level populations for specified ion.
        possible keyword arguments include temperature, eDensity, pDensity, radTemperature and rStar
        populate assumes that all of the population in the higher ionization stages exists only
        in the ground level
        use drPopulate() for cases where the population of various levels in the higher ionization stage
        figure into the calculation
        """
        nlvls = self.Nlvls

        nwgfa = self.Nwgfa
        nscups = self.Nscups
        npsplups = self.Npsplups
        nauto = self.Nauto
        # the following is the number of lower levels to ionize the current ion  from
        #  fixes the problem that the ionization rates are really just strictly for the ground leve
        # which can become depopulated with increased density


        if hasattr(self, 'PDensity'):
            protonDensity = self.PDensity
        else:
            self.p2eRatio()
            self.PDensity = self.ProtonDensityRatio*self.EDensity
            protonDensity = self.PDensity
            print(' proton density not specified, set to \"default\" ')
        #
        if hasattr(self, 'RadTemperature') and hasattr(self, 'RStar'):
            radTemperature = self.RadTemperature
            rStar = self.RStar
        rec = 0
        ci = 0
        if self.Nrrlvl or self.Nauto:
            rec = 1
        elif not popCorrect:
            rec = 0
        if verbose:
            print('rec = %5i'%(rec))
        if self.Nrrlvl:
            rrlvl = self.Rrlvl
            if hasattr(self, 'RrlvlRate'):
                rrlvlRate = self.RrlvlRate
            else:
                self.rrlvlDescale('rrlvl')
                rrlvlRate = self.RrlvlRate
            rrLvlIdx =[i for i,lvl1 in enumerate(self.Rrlvl['lvl1']) if lvl1 == 1]
            self.RrLvlIdx = rrLvlIdx
        temperature = self.Temperature
        eDensity = self.EDensity
        if rec:
            # get ionization rate of this current ion
            self.ionizRate()
            self.recombRate()
            #  get the higher ionization stage and its recombination rates to this ion
            highers = util.zion2name(self.Z, self.Ion+1)
            higher = ion(highers, temperature, abundance=self.Abundance, setup=0)
            higher.setupIonrec()
            higher.recombRate()
        #  the populating matrix for radiative transitions and autoionization transitions
        rad = np.zeros((nlvls+ci+rec,nlvls+ci+rec),np.float64)

        for iwgfa in range(nwgfa):
            l1 = self.Wgfa["lvl1"][iwgfa]-1
            l2 = self.Wgfa["lvl2"][iwgfa]-1
            rad[l1+ci,l2+ci] += self.Wgfa["avalue"][iwgfa]
            rad[l2+ci,l2+ci] -= self.Wgfa["avalue"][iwgfa]
            # photo-excitation and stimulated emission
            if self.RadTemperature:
                if not self.RStar:
                    dilute = 0.5
                else:
                    dilute = util.dilute(self.RStar)
                # next - don't include autoionization lines
                if abs(self.Wgfa['wvl'][iwgfa]) > 0.:
                    de = const.invCm2Erg*(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                    dekt = de/(const.boltzmann*self.RadTemperature)
                    # photoexcitation
                    phexFactor = dilute*(float(self.Elvlc['mult'][l2])/float(self.Elvlc['mult'][l1]))/(np.exp(dekt) -1.)
                    rad[l2+ci,l1+ci] += self.Wgfa["avalue"][iwgfa]*phexFactor
                    rad[l1+ci,l1+ci] -= self.Wgfa["avalue"][iwgfa]*phexFactor
                    # stimulated emission
                    stemFactor = dilute/(np.exp(-dekt) -1.)
                    rad[l1+ci,l2+ci] += self.Wgfa["avalue"][iwgfa]*stemFactor
                    rad[l2+ci,l2+ci] -= self.Wgfa["avalue"][iwgfa]*stemFactor

        # autoionization rates
        for iauto in range(nauto):
            l1 = self.Auto["lvl1"][iauto] - 1
            l2 = self.Auto["lvl2"][iauto] - 1
            # all autoionization eventually goes to the ground level of the higher ion
            rad[l1+ci+rec,l2+ci] += self.Auto["avalue"][iauto]
            rad[l2+ci,l2+ci] -= self.Auto["avalue"][iauto]


        if self.Nscups:
            self.upsilonDescale()
            #ups = self.Upsilon['upsilon']
            exRate = self.Upsilon['exRate']
            dexRate = self.Upsilon['dexRate']

        if self.Npsplups:
            self.upsilonDescale(prot=1)
            pexRate = self.PUpsilon['exRate']
            pdexRate = self.PUpsilon['dexRate']


        temp = self.Temperature
        ntemp = self.Ntemp
        eDensity = self.EDensity
        ndens = self.Ndens
        nscups = self.Nscups
#        cc = const.collision*self.EDensity
        #
        # (4 pi a0^2)^(3/2) = 6.6011e-24 (Badnell et al, 2003, A&A 406, 1151
        coef1 = 6.6011e-24*(const.hartree/(2.*const.boltzmann*self.Temperature))**1.5
        coef2 = (const.planck)**3/(2.*const.pi*const.emass*const.boltzmann*self.Temperature)**1.5

        #
#        if npsplups:
#            cp = const.collision*protonDensity
        #
        # get corrections for recombination and excitation
        #
        #
        # the way temperature and density are now (9/2015) handled as arrays of the same size
        # one the ndens == ntemp =1 case and the ndens >1 and ntemp>1 case are really needed
        #
        errorMessage = []
        # ---------------------------------------------------------------------------------
        #  --------------now for all values of ntemp and ndens -------------------------------
        # --------------------------------------------------------------------------------

        if verbose:
            print(' doing both ntemp: %5i  ndens:  %5i'%(ntemp, ndens))
        pop = np.zeros((self.NTempDens,nlvls),np.float64)
#            drPop = np.zeros((ntemp,nlvls),np.float64)
        fullPop = np.zeros((self.NTempDens, ci + nlvls + rec), np.float64)
        if self.Nrrlvl:
            rrTot = np.zeros(self.NTempDens, np.float64)
#            recRate = np.zeros((ntemp,nlvls),np.float64)
        #drEffRateTot = np.zeros(ntemp, np.float64)
        rrTot = np.zeros(self.NTempDens, np.float64)
        recTot = np.zeros(self.NTempDens, np.float64)
#
#
        higherPop = []
        for itemp in range(self.NTempDens):
            popmat = np.copy(rad)
            temp = self.Temperature[itemp]
            dens = self.EDensity[itemp]

            for iscups in range(nscups):
                l1 = self.Scups["lvl1"][iscups]-1
                l2 = self.Scups["lvl2"][iscups]-1
                popmat[l1+ci,l2+ci] += self.EDensity[itemp]*dexRate[iscups, itemp]
                popmat[l2+ci,l1+ci] += self.EDensity[itemp]*exRate[iscups, itemp]
                popmat[l1+ci,l1+ci] -= self.EDensity[itemp]*exRate[iscups, itemp]
                popmat[l2+ci,l2+ci] -= self.EDensity[itemp]*dexRate[iscups, itemp]

            # proton rates
            for ipslups in range(npsplups):
                l1 = self.Psplups["lvl1"][ipslups]-1
                l2 = self.Psplups["lvl2"][ipslups]-1
                popmat[l1+ci,l2+ci] += self.PDensity[itemp]*pdexRate[ipslups, itemp]
                popmat[l2+ci,l1+ci] += self.PDensity[itemp]*pexRate[ipslups, itemp]
                popmat[l1+ci,l1+ci] -= self.PDensity[itemp]*pexRate[ipslups, itemp]
                popmat[l2+ci,l2+ci] -= self.PDensity[itemp]*pdexRate[ipslups, itemp]

            if rec:
                for ilvl in range(0, self.GrndLevels):
                    popmat[-1,  ci + ilvl] += self.EDensity[itemp]*self.IonizRate['rate'][itemp]
                    popmat[ci + ilvl, ci + ilvl] -= self.EDensity[itemp]*self.IonizRate['rate'][itemp]

                if self.Nrrlvl:
#                    rrlvl = self.Rrlvl
#                    if hasattr(self, 'RrlvlRate'):
#                        rrlvl = self.Rrlvl
                    if not hasattr(self, 'RrlvlRate'):
                        self.rrlvlDescale()
                        # only include rr from ground level
                    for itrans in rrLvlIdx:
                        lvl1 = rrlvl['lvl1'][itrans]-1
                        lvl2 = rrlvl['lvl2'][itrans]-1
                        popmat[lvl2+ci, -1] += self.EDensity[itemp]*self.RrlvlRate['rate'][itrans, itemp]
                        popmat[-1, -1] -= self.EDensity[itemp]*self.RrlvlRate['rate'][itrans,itemp]
                        rrTot[itemp] += self.RrlvlRate['rate'][itrans, itemp]
                else:
                    rrTot[itemp] = 0.
                if self.Nauto:
                    autoLvl2 = []
                    for i, avalue in enumerate(self.Auto['avalue']):
                        l1 = self.Auto['lvl1'][i] - 1
                        l2 = self.Auto['lvl2'][i] - 1
                        autoLvl2.append(l2)
                        upperIdx = higher.Elvlc['lvl'].index(self.Auto['lvl1'][i])
                        gUpper = float(higher.Elvlc['mult'][l1])
                        gLower = float(self.Elvlc['mult'][l2])

                        ecm2 = self.Elvlc['ecm'][l2]
                        if ecm2 < 0.:
                            ecm2 = self.Elvlc['ecmth'][l2]
                        de1 = ecm2*const.invCm2Erg - self.Ip*const.ev2Erg
                        dekt1 = de1/(const.boltzmann*self.Temperature[itemp])
                        expkt = np.exp(-dekt1)

                        # only including dielectronic recombination from the lowest level
                        if l1 == 0:
                            dielRate = coef2[itemp]*gLower*expkt*avalue/(2.*gUpper)
                            popmat[ci + l2, -1] += self.EDensity[itemp]*dielRate
                            popmat[-1, -1] -= self.EDensity[itemp]*dielRate




                recTot[itemp] = rrTot[itemp]
                if recTot[itemp] < higher.RecombRate['rate'][itemp]:
                    popmat[ci, -1] += self.EDensity[itemp]*(higher.RecombRate['rate'][itemp] - recTot[itemp])
                    popmat[-1, -1] -= self.EDensity[itemp]*(higher.RecombRate['rate'][itemp] - recTot[itemp])
                if verbose:
                    print('itemp  rrTot         %5i %10.2e '%(itemp, rrTot[itemp]))
                    print('itemp, recTot RecombRate  %5i %10.2e %10.2e'%(itemp, recTot[itemp], higher.RecombRate['rate'][itemp]))
            else:
                recTot[itemp] = 0.

            # in this case, haven't completely accounted for recombination


#                    if verbose:
#                        print(' T recombRate, rrTot, rrTot %10.2e %10.2e %10.2e %10.2e'%(temp, higher.RrRate['rate'][itemp], rrTot[itemp], higher.DrRate['rate'][itemp]))

            norm = np.ones(nlvls+ci+rec,np.float64)
            self.popmat = np.copy(popmat)
            if ci:
                norm[0] = 0.
            if rec:
                norm[-1] = 0.
            popmat[nlvls+ci+rec-1] = norm
            b = np.zeros(nlvls+ci+rec,np.float64)
            b[nlvls+ci+rec-1] = 1.
            try:
                thispop = np.linalg.solve(popmat,b)
                if rec:
                    higherPop.append(thispop[ci+nlvls])
                else:
                    higherPop.append(-1)
                pop[itemp] = thispop[ci:ci+nlvls]
                fullPop[itemp] = thispop
            except np.linalg.LinAlgError:
                pop[itemp] = np.zeros(nlvls, np.float64)
                errorMessage.append('linealgError for T index %5i'%(itemp))
        #
            pop = np.where(pop > 0., pop, 0.)
        self.Population = {"temperature":temperature,"eDensity":eDensity,"population":pop,
            "protonDensity":protonDensity, "ci":ci, "rec":rec, 'popmat':popmat, 'fullPop':fullPop,
            'method':'populate',  'higherPop':higherPop, 'rrTot':rrTot}
        if len(errorMessage) > 0:
            self.Population['errorMessage'] = errorMessage


    def p2eRatio(self):
        """
        Calculates the proton density to electron density ratio using Eq. 7 of [1]_.

        Notes
        ------
        Uses the abundance and ionization equilibrium.

        References
        ----------
        .. [1] Young, P. R. et al., 2003, ApJS, `144, 135 <http://adsabs.harvard.edu/abs/2003ApJS..144..135Y>`_
        """
        if hasattr(self, 'Temperature'):
            temperature = self.Temperature
        else:
            temperature = self.IoneqAll['ioneqTemperature']
#        if not hasattr(self, 'AbundanceName'):
#            AbundanceName = self.Defaults['abundfile']
#        else:
#            AbundanceName = self.AbundanceName

        tmp_abundance = self.AbundAll
#        abundance = tmp_abundance['abundance'][tmp_abundance['abundance']>0]
        abundance = tmp_abundance[tmp_abundance > 0]
        denominator = np.zeros(len(self.IoneqAll['ioneqTemperature']))
        for i in range(len(abundance)):
            for z in range(1,i+2):
                denominator += z*self.IoneqAll['ioneqAll'][i,z,:]*abundance[i]

        p2eratio = abundance[0]*self.IoneqAll['ioneqAll'][0,1,:]/denominator
        nots = splrep(np.log10(self.IoneqAll['ioneqTemperature']),p2eratio,s=0)
        self.ProtonDensityRatio = splev(np.log10(temperature),nots,der=0,ext=1)



    def popPlot(self, top=10, levels=[], scale=0, plotFile=0, outFile=0, pub=0, addTitle=None,
        addLegend=True):
        """
        Plots populations vs temperature or eDensity.

        top specifies the number of the most highly populated levels to plot (the default)

        or can set levels to an array such as a list to set the desired levels to plot

        if scale is set, then the population, if plotted vs. density, is divided by density -
        only useful if plotting level populations vs density

        if pub is set, the want publication plots (bw, lw=2).

        if addLegend is set, a matplotlib legend is added
        """

        if pub:
            fontsize = 16
        else:
            fontsize = 14

        if hasattr(self, 'Population'):
            temperature = self.Population["temperature"]
            eDensity = self.Population["eDensity"]
            pop = self.Population["population"]
        else:
            self.populate()
            temperature = self.Population["temperature"]
            eDensity = self.Population["eDensity"]
            pop = self.Population["population"]

        #  for case of only a single density and temperature
        if len(pop.shape) == 1:
            spop = np.sort(pop)
            idx = np.argsort(pop)
            minPop = spop[-top:].min()/2.
            if top > pop.size:
                top = pop.size
            for itop in range(1, top+1):
                x = [idx[-itop], idx[-itop], idx[-itop]+1, idx[-itop]+1]
                y = [minPop, spop[-itop], spop[-itop], minPop]
                plt.semilogy(x, y, 'k')
            plt.axis([0, max(idx[-top:])+1, minPop, 1.])
            plt.xlabel('Level', fontsize=fontsize)
            plt.ylabel('Population', fontsize=fontsize)
            return

        # find the top most populated levels
        doTop = 0
        doLvl = 0
        if top and len(levels) >1:
            print(' to specifiy the levels to be plotted, top must be set to False or 0')
            return
        if top:
            doTop = 1
            lvl = self.Elvlc["lvl"]
            nlvls = self.Nlvls
            if top > nlvls:
                top = nlvls
            maxpop = np.zeros(nlvls,np.float64)
            for ilvl in range(nlvls):
                maxpop[ilvl] = pop[:,ilvl].max()
            lvlsort = np.take(lvl,np.argsort(maxpop))
            toplvl = lvlsort[-top:]
        elif len(levels) > 0:
            doLvl = 1
            toplvl = levels
            top = len(levels)
        ntemp = self.Ntemp
        ndens = self.Ndens

        ylabel = 'Population'
        title = self.Spectroscopic
        plt.figure()
        plt.ion()
        if ndens == 1:
            toppops = np.zeros((top, ntemp), np.float64)
            for ilvl in range(top):
                toppops[ilvl] = pop[:, toplvl[ilvl]-1]
            nonzero = toppops > 0.
            ymin = min(toppops[nonzero])
            for lvl in toplvl:
                # for some low temperature, populations can not be calculated
                good = pop[:, lvl-1] > 0
                if pub:
                    plt.loglog(temperature[good],pop[good,lvl-1], 'k',lw=2, label=str(lvl))
                else:
                    plt.loglog(temperature[good],pop[good,lvl-1], label=str(lvl))
                skip = 3
                if good.sum() == ntemp:
                    start = divmod(lvl,ntemp)[1]
                    for itemp in range(start,ntemp,ntemp//skip):
                        plt.text(temperature[itemp],pop[itemp,lvl-1],str(lvl))
                else:
                    newtemp = []
                    for i, one in enumerate(temperature):
                        if good[i]:
                            newtemp.append(one)
                    start = divmod(lvl, len(newtemp))[1] + ntemp - good.sum()
                    for itemp in range(start,ntemp,ntemp//skip):
                        plt.text(temperature[itemp],pop[itemp,lvl-1],str(lvl))
            xlabel = 'Temperature (K)'
            plt.xlabel(xlabel,fontsize=fontsize)
            plt.ylabel(ylabel,fontsize=fontsize)
            if addTitle is None:
                dstr = ' -  Density = %10.2e (cm$^{-3}$)' % eDensity[0]
            else:
                dstr = ' -  Density = %10.2e (cm$^{-3}$) %s' %(eDensity[0], addTitle)

            plt.title(title+dstr,fontsize=fontsize)
            plt.xlim(temperature.min(),temperature.max())
            if doTop:
                plt.ylim(ymin,1.2)
            elif doLvl:
                ymin = 2.
                ymax = 0.
                for ilvl in levels:
                    pretty = self.Elvlc['pretty'][ilvl-1]
                    print(' ilvl %5i %s min max %12.2e %12.2e'%(ilvl, pretty, pop[:, ilvl-1].min(), pop[:, ilvl-1].max()))
                    if pop[:, ilvl-1].max() > ymax:
                        ymax = pop[:, ilvl-1].max()
                    if pop[:, ilvl-1].min() < ymin:
                        ymin =  pop[:, ilvl-1].min()
                print(' doLvl ymin %12.2e  ymax %12.2e'%(ymin, ymax))
                plt.ylim(ymin/1.1, ymax*1.1 )
            if addLegend:
                plt.legend(loc = 'lower right')

        elif ntemp == 1:
            xlabel = r'Electron Density (cm$^{-3}$)'
            toppops = np.zeros((top, ndens), np.float64)
            for ilvl in range(top):
                toppops[ilvl] = pop[:, toplvl[ilvl]-1]
            nonzero = toppops > 0.
            ymin = min(toppops[nonzero])
            for lvl in toplvl:
                 # for some low temperature, populations can not be calculated
                good = pop[:, lvl-1] > 0.
                if scale:
                    plt.loglog(eDensity[good],pop[good,lvl-1]/eDensity[good], lw=1.5, label=str(lvl))
                elif pub:
                    plt.loglog(eDensity[good],pop[good,lvl-1], 'k', lw=2, label=str(lvl))
                else:
                    plt.loglog(eDensity[good],pop[good,lvl-1], label=str(lvl))
                skip = 3
                if good.sum() == ndens:
                    start = divmod(lvl,ndens)[1]
                    for idens in range(start,ndens,ndens//skip):
                        if scale:
                            plt.text(eDensity[idens],pop[idens,lvl-1]/eDensity[idens],str(lvl))
                        else:
                            plt.text(eDensity[idens],pop[idens,lvl-1],str(lvl))
                else:
                    newdens = []
                    for i, one in enumerate(eDensity):
                        if good[i]:
                            newdens.append(one)
                    start = divmod(lvl, len(newdens))[1] + ndens - good.sum()
                    for idens in range(start,ndens,ndens//skip):
                        plt.text(eDensity[idens],pop[idens, lvl-1],str(lvl))
            plt.xlabel(xlabel,fontsize=fontsize)
            plt.ylabel(ylabel,fontsize=fontsize)
            if addTitle is None:
                tstr = ' -  T = %10.2e (K)' % temperature[0]
            else:
                tstr = ' -  T = %10.2e (K) %s' % (temperature[0] , addTitle)
            plt.title(title+tstr,fontsize=fontsize)
            plt.xlim(eDensity[eDensity.nonzero()].min(),eDensity.max())
            if doTop:
                plt.ylim(ymin,1.2)
            elif doLvl:
                ymin = 2.
                ymax = 0.
                for ilvl in levels:
                    pretty = self.Elvlc['pretty'][ilvl-1]
                    if scale:
                        print(' ilvl %5i %s min max %12.2e %12.2e'%(ilvl, pretty, pop[:, ilvl-1].min(), pop[:, lvl-1].max()))
                        if (pop[:, ilvl-1]/eDensity).max() > ymax:
                            ymax = (pop[:, ilvl-1]/eDensity).max()
                        if (pop[:, ilvl-1]/eDensity).min() < ymin:
                            ymin =  (pop[:, ilvl-1]/eDensity).min()
                    else:
                        print(' ilvl %5i %s min max %12.2e %12.2e'%(ilvl, pretty, pop[:, ilvl-1].min(), pop[:, lvl-1].max()))
                        if pop[:, ilvl-1].max() > ymax:
                            ymax = pop[:, ilvl-1].max()
                        if pop[:, ilvl-1].min() < ymin:
                            ymin =  pop[:, ilvl-1].min()
                print(' doLvl ymin %12.2e  ymax %12.2e'%(ymin, ymax))
                plt.ylim(ymin/1.1, ymax*1.1 )
            if addLegend:
                plt.legend(loc='lower right')
        else:
            ax = plt.subplot(111)
            toppops = np.zeros((top, ntemp), np.float64)
            for ilvl in range(top):
                toppops[ilvl] = pop[:, toplvl[ilvl]-1]
            nonzero = toppops > 0.
            ymin = min(toppops[nonzero])
            for lvl in toplvl:
                # for some low temperature, populations can not be calculated
                good = pop[:, lvl-1] > 0
                if pub:
                    plt.loglog(temperature[good],pop[good,lvl-1], 'k', lw=2)
                else:
                    plt.loglog(temperature[good],pop[good,lvl-1])
                skip = 3
                if good.sum() == ntemp:
                    start = divmod(lvl,ntemp)[1]
                    for itemp in range(start,ntemp,ntemp//skip):
                        plt.text(temperature[itemp],pop[itemp,lvl-1],str(lvl))
                else:
                    newtemp = []
                    for i, one in enumerate(temperature):
                        if good[i]:
                            newtemp.append(one)
                    start = divmod(lvl, len(newtemp))[1] + ntemp - good.sum()
                    for itemp in range(start,ntemp,ntemp//skip):
                        plt.text(temperature[itemp],pop[itemp,lvl-1],str(lvl))
            xlabel = 'Temperature (K)'
            plt.xlabel(xlabel,fontsize=fontsize)
            plt.ylabel(ylabel,fontsize=fontsize)
            plt.axis([temperature.min(),temperature.max(), ymin, 1.2])
            plt.text(0.1, 0.5,title, horizontalalignment='center', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            #
            ax2 = plt.twiny()
            xlabel = r'Electron Density (cm$^{-3}$)'
            plt.xlabel(xlabel, fontsize=fontsize)
            plt.loglog(eDensity,pop[:,toplvl[0]], visible=False)
            ax2.xaxis.tick_top()
        plt.tight_layout()
        if outFile:
            plt.savefig(outFile)
        self.Population['toplvl'] = toplvl

    def emiss(self, allLines=True):
        """
        Calculate the emissivities for lines of the specified ion.

        units:  ergs s^-1 str^-1

        Does not include elemental abundance or ionization fraction

        Wavelengths are sorted

        set allLines = True to include unidentified lines
        """

        if hasattr(self, 'Population'):
            pop = self.Population['population']
        else:
            self.populate()
            pop = self.Population["population"]

        wvl = np.asarray(self.Wgfa["wvl"], np.float64)
        obs = np.where(wvl > 0., 'Y', 'N')

#        if self.Defaults['wavelength'] == 'ev':
#            wvl = const.ev2Ang/np.asarray(wvl)
#        if self.Defaults['wavelength'] == 'kev':
#            wvl = 1.e-3*const.ev2Ang/np.asarray(wvl)
#        elif self.Defaults['wavelength'] == 'nm':
#            wvl = wvl/10.

        if allLines:
            wvl=np.abs(wvl)
        l1  =  np.asarray(self.Wgfa['lvl1'], 'int64')
        l2 = np.asarray(self.Wgfa["lvl2"], 'int64')
        avalue = np.asarray(self.Wgfa["avalue"], np.float64)
        if 'pretty1' in self.Wgfa.keys():
            pretty1 = np.asarray(self.Wgfa['pretty1'])
        if 'pretty2' in self.Wgfa.keys():
            pretty2 = np.asarray(self.Wgfa['pretty2'])

        # make sure there are lines in the wavelength range, if specified

        # two-photon decays have wvl=0 and nonzero avalues
        nonzed = wvl != 0.
        wvl = wvl[nonzed]
        l1 = l1[nonzed]
        l2 = l2[nonzed]
        avalue = avalue[nonzed]
        pretty1 = pretty1[nonzed]
        pretty2 = pretty2[nonzed]
        obs = obs[nonzed]
        nwvl = len(wvl)

#        wvl = np.asarray(self.Wgfa["wvl"], np.float64)
#        obs = np.where(wvl > 0., 'Y', 'N')

        if self.Defaults['wavelength'] == 'ev':
            wvl = const.ev2Ang/np.asarray(wvl)
        if self.Defaults['wavelength'] == 'kev':
            wvl = 1.e-3*const.ev2Ang/np.asarray(wvl)
        elif self.Defaults['wavelength'] == 'nm':
            wvl = wvl/10.

#        if allLines:
#            wvl=np.abs(wvl)
#        l1  =  np.asarray(self.Wgfa['lvl1'], 'int64')
#        l2 = np.asarray(self.Wgfa["lvl2"], 'int64')
#        avalue = np.asarray(self.Wgfa["avalue"], np.float64)
#        if 'pretty1' in self.Wgfa.keys():
#            pretty1 = np.asarray(self.Wgfa['pretty1'])
#        if 'pretty2' in self.Wgfa.keys():
#            pretty2 = np.asarray(self.Wgfa['pretty2'])

        try:
            ntempden,nlvls = pop.shape
            em = np.zeros((nwvl, ntempden),np.float64)
        except:
            ntempden = 1
            em = np.zeros(nwvl,np.float64)

#        plotLabels = {}
#        if self.Defaults['wavelength'] == 'angstrom':
#            plotLabels["xLabel"] = "Angstroms"
#        elif self.Defaults['wavelength'] == 'nm':
#            plotLabels["xLabel"] = "nanometers"
#        elif self.Defaults['wavelength'] == 'kev':
#            plotLabels["xLabel"] = "kev"
        xlabel = self.Labels['xlabel']

        if self.Defaults['flux'] == 'energy':
            factor = const.planck*const.light/(4.*const.pi*1.e-8*wvl)
#            plotLabels["yLabel"] = "ergs cm^-3 s^-1"
        elif self.Defaults['flux'] == 'photon':
            factor = np.ones((nwvl),np.float64)/(4.*const.pi)
#            plotLabels["yLabel"] = "photons cm^-3 s^-1"
        ylabel = self.Labels['emissYlabel']
#        if ntempden > 1:
        for itempden in range(ntempden):
            for iwvl in range(nwvl):
                p = pop[itempden,l2[iwvl]-1]
                em[iwvl, itempden] = factor[iwvl]*p*avalue[iwvl]
#        else:
#            for iwvl in range(0,nwvl):
#                p = pop[l2[iwvl]-1]
#                em[iwvl] = factor[iwvl]*p*avalue[iwvl]
#            if self.Defaults['wavelength'] == 'kev':
#                wvl = const.ev2Ang/np.asarray(wvl)
#            elif self.Defaults['wavelength'] == 'nm':
#                wvl = wvl/10.
        nlvl = len(l1)
        ionS = np.asarray([self.IonStr]*nlvl)
        Emiss = {'ionS':ionS,"wvl":wvl, "emiss":em, 'lvl1':l1, 'lvl2':l2, 'avalue':avalue,
            'obs':obs, 'pretty1':pretty1, 'pretty2':pretty2, 'xlabel':xlabel,  'ylabel':ylabel}
        self.Emiss = Emiss
        return

    def emissList(self, index=None, wvlRange=None, wvlRanges=None, top=10, relative=0, outFile=0 ):
        """
        List the emissivities.

        wvlRange, a 2 element tuple, list or array determines the wavelength range

        Top specifies to plot only the top strongest lines, default = 10

        normalize = 1 specifies whether to normalize to strongest line, default = 0
        """
        #
        if outFile:
            output = open(outFile, 'w')
        #
        if not hasattr(self, 'Emiss'):
            try:
                self.emiss()
            except:
                print(' emissivities not calculated or emiss() is unable to calculate them')
                print(' perhaps the temperature and/or eDensity are not set')
                return
        #
        #
        temperature = self.Temperature
        eDensity = self.EDensity
        #
        ndens = self.Ndens
        ntemp = self.Ntemp
        #
        if ndens == 1 and ntemp == 1:
            dstr = ' -  Density = %10.2e (cm$^{-3}$)' %(eDensity)
            tstr = ' -  T = %10.2e (K)' %(temperature)
            index = 0
        elif ndens == 1 and ntemp > 1:
            if index is None:
                index = ntemp//2
            print('using index = %5i specifying temperature =  %10.2e'%(index, temperature[index]))
            self.Message = 'using index = %5i specifying temperature  =   %10.2e'%(index, temperature[index])
        elif ndens > 1 and ntemp == 1:
            if index is None:
                index = ntemp//2
            print('using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index]))
            self.Message = 'using index = %5i specifying eDensity = %10.2e'%(index, eDensity[index])
        elif ndens > 1 and ntemp > 1:
            if index is None:
                index = ntemp//2
            print('using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index]))
            self.Message = 'using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index])

        emissivity = copy.copy(self.Emiss)
        emiss = emissivity['emiss'][:, index]
        ionS = emissivity['ionS']
        wvl = emissivity['wvl']
        lvl1 = emissivity['lvl1']
        lvl2 = emissivity['lvl2']
        avalue = emissivity['avalue']
        obs = emissivity['obs']
        pretty1 = emissivity['pretty1']
        pretty2 = emissivity['pretty2']
#        plotLabels = emissivity['plotLabels']
        #
        if wvlRange:
            wvlIndex = util.between(wvl,wvlRange)
        elif wvlRanges:
            wvlIndex = []
            for awvlRange in wvlRanges:
                wvlIndex.extend(util.between(wvl,awvlRange))
        else:
            wvlIndex = range(wvl.size)
        #
        #
        # get those in the right wavelength range
        #
        emiss = emiss[wvlIndex]
        ionS = ionS[wvlIndex]
        wvl = wvl[wvlIndex]
        lvl1 = lvl1[wvlIndex]
        lvl2 = lvl2[wvlIndex]
        avalue = avalue[wvlIndex]
        obs = obs[wvlIndex]
        pretty1 = pretty1[wvlIndex]
        pretty2 = pretty2[wvlIndex]
        #
        self.Error = 0
        if wvl.size == 0:
            print('No lines in this wavelength interval')
            self.Error = 1
            self.Message = 'No lines in this wavelength interval'
            return
        #
        elif top == 0:
            top = wvl.size
        elif top > wvl.size:
            top = wvl.size
        #
        # now sort by intensity/emissivity
        isrt = np.argsort(emiss)
        ionS = ionS[isrt[-top:]]
        wvl = wvl[isrt[-top:]]
        lvl1 = lvl1[isrt[-top:]]
        lvl2 = lvl2[isrt[-top:]]
        obs = obs[isrt[-top:]]
        emiss = emiss[isrt[-top:]]
        avalue = avalue[isrt[-top:]]
        pretty1 = pretty1[isrt[-top:]]
        pretty2 = pretty2[isrt[-top:]]
        #
    # must follow setting top
        #
        if relative:
            emiss = emiss/emiss[:top].max()
        #
        #
        idx = np.argsort(wvl)
        #
        fmt = '%5s %5i %5i %25s - %25s %12.3f %12.3e %12.2e %1s'
        fmtTitle = '%5s %5s %5s %25s - %25s %12s %12s %12s %1s'
        print( '  ')
        print( '------------------------------------------')
        print('  ')
#        print(' Ion   lvl1  lvl2         lower                     upper                   Wvl(A)   Emissivity      A value Obs')
        title = fmtTitle%('Ion', 'lvl1', 'lvl2', 'lower', 'upper', self.Labels['listXlabel'],
            'Emissivity', 'A value', 'Obs')
        print(title)
        for kdx in idx:
            print(fmt%(ionS[kdx], lvl1[kdx], lvl2[kdx], pretty1[kdx], pretty2[kdx], wvl[kdx],
                emiss[kdx], avalue[kdx], obs[kdx]))
        print('   ')
        print(' ------------------------------------------')
        print('   ')
        #
#        self.Emiss['wvlTop'] = wvl[idx]
#        self.Emiss['emissTop'] = emiss[idx]
#        self.Emiss['idx'] = idx
        if outFile:
            fmt = '%5s %5i %5i %25s - %25s %12.3f %12.3e %12.2e %1s'
            output.write('   \n')
            output.write(' ------------------------------------------ \n')
            output.write('   ')
            output.write(' Ion   lvl1  lvl2         lower                     upper                   Wvl(A)   Emissivity      A value Obs  \n')
            for kdx in idx:
                output.write(fmt%(ionS[kdx], lvl1[kdx], lvl2[kdx], pretty1[kdx], pretty2[kdx], wvl[kdx], emiss[kdx], avalue[kdx], obs[kdx]))
            output.write('    \n')
            output.write(' ------------------------------------------  \n')
            output.write('   \n ')
            output.close()
        #
        # ---------------------------------------------------------------------------
        #
    def emissPlot(self, index=None, wvlRange=None,  top=10, linLog='lin', relative=0,  verbose=0,
        plotFile = 0, saveFile=0 ):
        '''Plot the emissivities.

        wvlRange, a 2 element tuple, list or array determines the wavelength range

        Top specifies to plot only the top strongest lines, default = 10

        linLog specifies a linear or log plot, want either lin or log, default = lin

        normalize = 1 specifies whether to normalize to strongest line, default = 0'''
        #
        fs = 14
        title = self.Spectroscopic
        #
        if hasattr(self, 'Emiss'):
            em = self.Emiss
        else:
            try:
                self.emiss()
                em = self.Emiss
            except:
                print(' emissivities not calculated and emiss() is unable to calculate them')
                print(' perhaps the temperature and/or eDensity are not set')
                return
        emiss = em['emiss']
        wvl = em['wvl']
        temperature = self.Temperature
        eDensity = self.EDensity
        #
        ndens = self.Ndens
        ntemp = self.Ntemp
        #
        if ndens == 1 and ntemp == 1:
            dstr = ' -  Density = %10.2e (cm$^{-3}$)' %(eDensity)
            tstr = ' -  T = %10.2e (K)' %(temperature)
            index = 0
        elif ndens == 1 and ntemp > 1:
            if index is None:
                index = ntemp//2
                print('using index = %5i specifying temperature =  %10.2e'%(index, temperature[index]))
                self.Message = 'using index = %5i specifying temperature = %10.2e'%(index, temperature[index])
            dstr = ' -  Density = %10.2e (cm$^{-3}$)' % eDensity[index]
            tstr = ' -  T = %10.2e (K)' % temperature[index]
#            emiss = emiss[:, index]
        elif ndens > 1 and ntemp == 1:
            if index is None:
                index = ntemp//2
                print('using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index]))
                self.Message = 'using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index])
            emiss = emiss[:, index]
            dstr = ' -  Density = %10.2e (cm$^{-3}$)' % eDensity[index]
            tstr = ' -  T = %10.2e (K)' % temperature
#            emiss = emiss[:, index]
        elif ndens > 1 and ntemp > 1:
            if index is None:
                index = ntemp//2
                print('using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index]))
                self.Message = 'using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index])
#            emiss = emiss[:, index]
            dstr = ' -  Density  =  %10.2e (cm$^{-3}$)' % eDensity[index]
            tstr=' -  T = %10.2e (K)' % temperature[index]

        emiss = emiss[:, index]

        if wvlRange is not None:
            wvlIndex = util.between(wvl, wvlRange)
        else:
            wvlIndex = range(wvl.size)
        emiss = emiss[wvlIndex]
        wvl = wvl[wvlIndex]
        #
        self.Error = 0
        if wvl.size == 0:
            print('No lines in this wavelength interval')
            self.Error = 1
            self.Message = 'No lines in this wavelength interval'
            return
        elif top == 0:
            top = wvl.size
        elif wvl.size > top:
            isrt = np.argsort(emiss)
            wvl = wvl[isrt[-top:]]
            emiss = emiss[isrt[-top:]]
        else:
            top = wvl.size
        # must follow setting top
        #
        plt.figure()
        ylabel = self.Labels['emissYlabel']
        if relative:
            emiss = emiss/emiss[:top].max()
            ylabel += ' (Relative)'
        #
#        xlabel = 'Wavelength ('+self.Defaults['wavelength'] +')'
        xlabel = self.Labels['xlabel']
        #
        ymin = 10.**(np.log10(emiss.min()).round(0)-0.5 )
        #
        plt.ion()
        #
        for idx in range(top):
            xx = [wvl[idx], wvl[idx]]
            if linLog == 'lin':
                yy = [0., emiss[idx]]
                plt.plot(xx, yy)
            else:
                yy = [ymin/10., emiss[idx]]
                plt.semilogy(xx, yy)
        plt.xlabel(xlabel, fontsize=fs)
        plt.ylabel(ylabel, fontsize=fs)
        plt.title(title+tstr+dstr, fontsize=fs)
        plt.tight_layout()
        if wvlRange:
            plt.axis([wvlRange[0], wvlRange[1], 0., 1.1*emiss.max()])
        if plotFile:
            plt.savefig(plotFile)
        #
        # -------------------------------------------------------------------------------------
        #
    def emissRatio(self, wvlRange=None, wvlRanges=None, top=10):
        """
        Plot the ratio of 2 lines or sums of lines.
        Shown as a function of density and/or temperature.
        For a single wavelength range, set wvlRange = [wMin, wMax]
        For multiple wavelength ranges, set wvlRanges = [[wMin1,wMax1],[wMin2,wMax2], ...]
        A plot of relative emissivities is shown and then a dialog appears for the user to
        choose a set of lines.
        """

        if hasattr(self, 'Emiss'):
            doEmiss = False
            em = self.Emiss
        else:
            doEmiss = True
        if doEmiss:
            # new values of temperature or eDensity
            self.emiss()
            em = self.Emiss
        fontsize = 14
        eDensity = self.EDensity
        temperature = self.Temperature
        ntemp = temperature.size

        if ntemp > 0:
            if temperature[0] == temperature[-1]:
                ntemp = 1
        ndens = eDensity.size
        if ndens > 0:
            if eDensity[0] == eDensity[-1]:
                ndens = 1
        print(' ndens = %5i ntemp = %5i'%(ndens, ntemp))

        emiss = em['emiss']
        ionS = em['ionS']
        wvl = em["wvl"]
        lineLabel = []
        for iline,  ions in enumerate(ionS):
            lineLabel.append(ions+' '+str(wvl[iline]))
        xlabel = self.Labels["xlabel"]
        ylabel = self.Labels["emissYlabel"]

        # find which lines are in the wavelength range if it is set
        if wvlRange:
            igvl = util.between(wvl,wvlRange)
        elif wvlRanges:
            igvl = []
            for awvlRange in wvlRanges:
                igvl.extend(util.between(wvl,awvlRange))
        else:
            igvl = range(len(wvl))
        #
        nlines = len(igvl)
        igvl = np.take(igvl,wvl[igvl].argsort())
        # find the top most intense lines
        if top > nlines:
            top = nlines
            #
        maxEmiss = np.zeros(nlines,np.float64)
        for iline in range(nlines):
            maxEmiss[iline] = emiss[igvl[iline]].max()
        for iline in range(nlines):
            if maxEmiss[iline] == maxEmiss.max():
                maxAll = emiss[igvl[iline]]

        igvlsort = np.take(igvl,np.argsort(maxEmiss))
        topLines = igvlsort[-top:]
        maxWvl = '%5.3f' % wvl[topLines[-1]]
        topLines = topLines[wvl[topLines].argsort()]

        # need to make sure there are no negative values before plotting
        good  =  np.where(emiss > 0.)
        emissMin = emiss[good].min()
        bad = np.where(emiss <= 0.)
        emiss[bad] = emissMin
        ntemp = self.Ntemp
        ndens = self.Ndens
        ylabel = 'Emissivity relative to '+maxWvl
        title = self.Spectroscopic

        if ndens==1 and ntemp==1:
            print(' only a single temperature and eDensity')
            return
        elif ndens == 1:
            xlabel = 'Temperature (K)'
            xvalues = self.Temperature
            outTemperature = self.Temperature
            outDensity = np.zeros(ntemp,np.float64)
            outDensity.fill(self.EDensity[0])
            desc_str = ' at  Density = %10.2e (cm)$^{-3}$' % self.EDensity[0]
        elif ntemp == 1:
            xvalues = self.EDensity
            outTemperature = np.zeros(ndens,np.float64)
            outTemperature.fill(self.Temperature[0])
            outDensity = self.EDensity
            xlabel = r'$\rm{Electron Density (cm)^{-3}}$'
            desc_str = ' at Temp = %10.2e (K)' % self.Temperature
        else:
            outTemperature = self.Temperature
            outDensity = self.EDensity
            xlabel = 'Temperature (K)'
            xvalues = self.Temperature
            desc_str = ' for variable Density'

        # put all actual plotting here
        plt.ion()
        #  maxAll is an array
        ymax = np.max(emiss[topLines[0]]/maxAll)
        ymin = ymax
        plt.figure()
        ax = plt.subplot(111)
        nxvalues = len(xvalues)

        # reversing is necessary - otherwise, get a ymin=ymax and a matplotlib error
        if self.Defaults['wavelength'] == 'angstrom' or self.Defaults['wavelength'] == 'nm':
            param0 = top-1
            param1 = -1
            param2 = -1
        elif self.Defaults['wavelength'] == 'ev' or self.Defaults['wavelength'] == 'kev':
            param0 = 0
            param1 = top
            param2 = 1

#        for iline in range(top-1, -1, -1):
#        for iline in range(top):
        for iline in range(param0, param1, param2):

            tline = topLines[iline]
            plt.loglog(xvalues,emiss[tline]/maxAll)
            if np.min(emiss[tline]/maxAll) < ymin:
                ymin = np.min(emiss[tline]/maxAll)
            if np.max(emiss[tline]/maxAll) > ymax:
                ymax = np.max(emiss[tline]/maxAll)
            skip = 2
            start = divmod(iline,nxvalues)[1]

            if self.Defaults['wavelength'] == 'angstrom':
                alabel = '%10.3f'%(wvl[tline])
            elif self.Defaults['wavelength'] == 'nm':
                alabel = '%10.4f'%(wvl[tline])
            elif self.Defaults['wavelength'] == 'ev':
                alabel = '%10.4e'%(wvl[tline])
            elif self.Defaults['wavelength'] == 'kev':
                alabel = '%10.3f'%(wvl[tline])


            for ixvalue in range(start,nxvalues,nxvalues//skip):
#                plt.text(xvalues[ixvalue],emiss[tline,ixvalue]/maxAll[ixvalue],str(wvl[tline]))
                plt.text(xvalues[ixvalue],emiss[tline,ixvalue]/maxAll[ixvalue],alabel)
        plt.xlim(xvalues.min(),xvalues.max())
        plt.xlabel(xlabel,fontsize=fontsize)
        plt.ylabel(ylabel,fontsize=fontsize)
        if ndens == ntemp and ntemp > 1:
            plt.text(0.07, 0.5,title, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            #
            ax2 = plt.twiny()
            xlabelDen = r'Electron Density (cm$^{-3}$)'
            plt.xlabel(xlabelDen, fontsize=fontsize)
            plt.loglog(eDensity,emiss[topLines[top-1]]/maxAll, visible=False)
            ax2.xaxis.tick_top()
            plt.ylim(ymin/1.2, 1.2*ymax)
        else:
            plt.ylim(ymin/1.2, 1.2*ymax)
            plt.title(title+desc_str,fontsize=fontsize)
        plt.tight_layout()
        plt.draw()
        #  need time to let matplotlib finish plotting
        time.sleep(0.5)
        #
        # get line selection  ************************************************************
        #
        selectTags = []
        for itop in topLines:

            if self.Defaults['wavelength'] == 'angstrom':
                alabel = '%s %10.3f'%(ionS[itop],  wvl[itop])
            elif self.Defaults['wavelength'] == 'nm':
                alabel = '%s %10.4f'%(ionS[itop],  wvl[itop])
            elif self.Defaults['wavelength'] == 'ev':
                alabel = '%s %10.4e'%(ionS[itop],  wvl[itop])
            elif self.Defaults['wavelength'] == 'kev':
                alabel = '%s %10.3f'%(ionS[itop],  wvl[itop])

#            selectTags.append(ionS[itop]+ ' '+ str(wvl[itop]))
            selectTags.append(alabel)
        #
        numden = chGui.gui.choice2Dialog(selectTags)
#        numden = gui.choice2Dialog(wvl[topLines])
        #
        # num_idx and den_idx are tuples
        #
        num_idx = numden.numIndex
        if len(num_idx) == 0:
            print(' no numerator lines were selected')
            return
        #
        den_idx = numden.denIndex
        if len(den_idx) == 0:
            print(' no denominator lines were selected')
            return
        #
        numEmiss = np.zeros(len(xvalues),np.float64)
        for aline in num_idx:
            numEmiss += emiss[topLines[aline]]
        #
        denEmiss = np.zeros(len(xvalues),np.float64)
        for aline in den_idx:
            denEmiss += emiss[topLines[aline]]
        #
        # plot the desired ratio
        #  maxAll is an array
        plt.figure()
        ax = plt.subplot(111)
        plt.loglog(xvalues,numEmiss/denEmiss)
        plt.xlim(xvalues.min(),xvalues.max())
        plt.xlabel(xlabel,fontsize=fontsize)
        plt.ylabel('Emissivity Ratio ('+self.Defaults['flux']+')',fontsize=fontsize)
        desc = ''
        for aline in num_idx:
            desc += ' ' + selectTags[aline]
#            desc += ' ' + str(wvl[topLines[aline]])
        desc += ' / '
        for aline in den_idx:
            desc += ' ' + selectTags[aline]
#            desc += ' ' + str(wvl[topLines[aline]])
        if ndens == ntemp and ntemp > 1:
            plt.text(0.07, 0.5,desc, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            #
            ax2 = plt.twiny()
            xlabelDen = r'Electron Density (cm$^{-3}$)'
            plt.xlabel(xlabelDen, fontsize=fontsize)
            plt.loglog(eDensity,numEmiss/denEmiss, visible=False)
            ax2.xaxis.tick_top()
        else:
#            plt.ylim(ymin, ymax)
            plt.title(desc,fontsize=fontsize)
        plt.tight_layout()
        #
        intensityRatioFileName = self.IonStr
        for aline in num_idx:
            intensityRatioFileName += '_%3i'%(wvl[topLines[aline]])
        intensityRatioFileName +='_2'
        for aline in den_idx:
            intensityRatioFileName += '_%3i'%(wvl[topLines[aline]])
        intensityRatioFileName += '.rat'
        self.IntensityRatio = {'ratio':numEmiss/denEmiss,'desc':desc,
                'temperature':outTemperature,'eDensity':outDensity,'filename':intensityRatioFileName, 'numIdx':num_idx, 'denIdx':den_idx}

    def intensity(self, allLines=1, verbose=0):
        """
        Calculate  the intensities for lines of the specified ion.

        units:  ergs cm^-3 s^-1 str^-1

        includes elemental abundance and ionization fraction.

        the emission measure 'em' is included if specified
        """

        # so we know that it has been applied
        if not hasattr(self, 'Emiss'):
            self.emiss(allLines=allLines)
            emiss = copy.copy(self.Emiss)
        else:
            emiss = copy.copy(self.Emiss)
        if 'errorMessage'  in emiss.keys():
            self.Intensity = {'errorMessage': emiss['errorMessage']}
            return

        # everything in emiss should be a numpy array
        emissivity = emiss['emiss']
        ionS = emiss['ionS']
        wvl = emiss['wvl']
        lvl1 = emiss['lvl1']
        lvl2 = emiss['lvl2']
        obs = emiss['obs']
        pretty1 = emiss['pretty1']
        pretty2 = emiss['pretty2']
        avalue = emiss['avalue']
        errorMessage = None

        xlabel = self.Labels['xlabel']
        ylabel = self.Labels['intensityYlabel']

        if np.array_equal(self.Em, np.ones_like(self.Em)):
            ylabel += '($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'

        if hasattr(self, 'Abundance'):
            ab = self.Abundance
        else:
            self.Abundance = io.abundanceRead()
            ab = self.Abundance

        if len(emissivity.shape) > 1:
            nwvl, ntempden =  emissivity.shape
            if verbose:
                print('nwvl:  %5i  ntempden:  %5i'%(nwvl, ntempden))
            intensity = np.zeros((ntempden, nwvl),np.float64)
            good = self.IoneqOne > 0.
            if good.sum() == 0.:
                intensity = np.zeros((ntempden, nwvl),np.float64)
                errorMessage = 'ioneq = zero in temperature range'
                self.Intensity = {'intensity':intensity, 'errorMessage':errorMessage}
                return
            else:
                thisIoneq = self.IoneqOne
                for it in range(ntempden):
                    intensity[it] = ab*thisIoneq[it]*emissivity[:, it]*self.Em[it]/self.EDensity[it]
        else:
            thisIoneq = self.IoneqOne
            nwvl = len(emissivity)
            ntempden = 1
            intensity = ab*thisIoneq*emissivity*self.Em/self.EDensity
        if ntempden == 1:
            integrated = intensity
        else:
            integrated = intensity.sum(axis=0)
        Intensity = {'intensity':intensity, 'integrated':integrated,'ionS':ionS, 'wvl':wvl, 'lvl1':lvl1,
            'lvl2':lvl2, 'pretty1':pretty1, 'pretty2':pretty2,  'obs':obs, 'avalue':avalue,
            'em':self.Em, 'xlabel':xlabel, 'ylabel':ylabel}
        if errorMessage is not None:
            Intensity['errorMessage'] = errorMessage
        self.Intensity = Intensity

    def boundBoundLoss(self, allLines=1):
        """
        Calculate  the summed radiative loss rate for all spectral lines of the specified ion.

        Parameters
        ----------

        allLines : `bool`
            If True, include losses from both observed and unobserved lines.
            If False, only include losses from observed lines.

        includes elemental abundance and ionization fraction.


        Returns
        -------

        creates the attribute:

        BoundBoundLoss : `dict` with the keys below.

            *rate* : the radiative loss rate (:math:`\mathrm{erg \, cm^{-3}} \, \mathrm{s}^{-1}`) per unit emission measure.

            *temperature* : (K).

            *eDensity* : electron density (:math:`\mathrm{cm^{-3}}`)


        """

        self.emiss(allLines=allLines)
        emiss = self.Emiss
        if 'errorMessage'  in emiss.keys():
            self.Intensity = {'errorMessage': self.Spectroscopic+' no lines in this wavelength region'}
            return
        em = emiss['emiss']
        wvl = emiss['wvl']
        eDensity = self.EDensity
        if hasattr(self, 'Abundance'):
            ab = self.Abundance
        else:
            self.Abundance = io.abundanceRead()
            ab = self.Abundance
        if hasattr(self, 'IoneqOne'):
            thisIoneq = self.IoneqOne
        else:
            self.ioneqOne()
            thisIoneq = self.IoneqOne
        # should probably replace this with an if statement based on len of em.shape
        try:
            nwvl, ntempden = em.shape
            intensity = np.zeros((ntempden, nwvl),np.float64)
            if thisIoneq.size == 1:
                thisIoneq = np.ones(ntempden, np.float64)*thisIoneq
            for it in range(ntempden):
                if self.Defaults['flux'] != 'energy':
                    intensity[it] = 4.*const.pi*(const.planck*const.light*1.e+8/wvl)*ab*thisIoneq[it]*em[:, it]
                else:
                    intensity[it] = 4.*const.pi*ab*thisIoneq[it]*em[:, it]/eDensity[it]
            loss = intensity.sum(axis=1)
        except:
            nwvl = len(em)
            ntempden = 1
            if self.Defaults['flux'] != 'energy':
                intensity = 4.*const.pi*(const.planck*const.light*1.e+8/wvl)*ab*thisIoneq*em
            else:
                intensity = 4.*const.pi*ab*thisIoneq*em/eDensity[it]
            loss = intensity.sum()
        self.BoundBoundLoss = {'rate':loss, 'temperature':self.Temperature, 'eDensity':self.EDensity}

    def intensityRatioInterpolate(self, data, scale='lin', plot=False, verbose=False):
        '''
        Take a set of data and interpolate against the IntensityRatio.

        Given a set of intensities ratio data it correlates them to electron density
        values. This is done by interpolating against the curve computed from
        the IntenstyRatio (dict). For the interpolation, the B-spline representation
        is used from scipy.

        Parameters
        ----------
        data : float or array-like
            The intensities ratio data to be converted into density values.

        scale : {'lin', 'loglog', 'logx', 'logy'}, str, default : 'lin'
            The scale to be used for the interpolation (default 'lin' indicates linear).

        plot : bool, default : False
            If True the input data will be plotted on top of the IntensityRatio
            curve.

        verbose : bool, default : False
            If True displays a message about the physical quantity that was used
            (i.e. Temperature or Density).
            Additionally, prints the interpolated data in pairs of data,value
            (i.e. input intensities ratio data and corresponding density values).

        Returns
        -------

        Creates the attribute:

        IntensityRatioInterpolated : dict
            Dictionary that contains the input data and the derived density values,
            with keys

            *data* : input intensities ratio (no units)

            *value* : electron density (:math:`\mathrm{cm^{-3}}`)


        Examples
        --------
        >>> import ChiantiPy.core as ch
        >>> temp = 10**5.9
        >>> dens = 10**(np.linspace(8,12,num=50,endpoint=True))
        >>> fe12 = ch.ion('fe_12',temperature=temp,eDensity=dens)
        >>> fe12.intensityRatio(wvlRange=[186.,196.])
        >>> # Select for example lines 186.887 and 195.119 from the GUI
        >>> int_ratios = np.linspace(0.1,0.8,20)
        >>> fe12.intensityRatioInterpolate(int_ratios)
        >>> dens_values = fe12.IntensityRatioInterpolated['value']

        >>> # Same setup but using a 2D array of input data
        >>> x = 10
        >>> y = 30
        >>> int_ratios_map = np.random.random((x,y))
        >>> fe12.intensityRatioInterpolate(int_ratios_map,verbose=True)
        >>> dens_values_map = fe12.IntensityRatioInterpolated['value'].reshape((x,y))

        '''

        # First, what variable to use
        if self.IntensityRatio['temperature'].max() > self.IntensityRatio['temperature'].min():
            x = self.IntensityRatio['ratio']
            y = self.IntensityRatio['temperature']
            if verbose:
                print('Using temperature with %i5 values'%(len(x)))
                print(' number of values')
        else:
            x = self.IntensityRatio['ratio']
            y = self.IntensityRatio['eDensity']
            if verbose:
                print('Using density with %i5 values'%(len(x)))
                print(' number of values')

        #
        if x[0] > x[-1]:
            x = sorted(x)
            sy = []
            for idx in range(len(y) -1, -1, -1):
                sy.append(y[idx])
        else:
            sy = y

        # If necessary, flatten the input data in order to handle N-D arrays
        if len(np.shape(data)) > 1:
            if verbose:
                print('Input data have dimensions:')
                print(np.shape(data))
                print('They will be flattened!\n')
            data = np.array(data).flatten() # In case the data type is list

        #
        if 'lin' in scale:
            y2 = splrep(x, sy, s=0)
            interpolatedData = splev(data,y2)
            if plot:
                plt.plot(sy, x)
                plt.plot(interpolatedData, data, 'bD')
        elif scale == 'loglog':
            y2 = splrep(np.log(x), np.log(sy), s=0)
            interpolatedData = np.exp(splev(np.log(data),y2))
            if plot:
                plt.loglog(sy, x)
                plt.loglog(interpolatedData, data, 'bD')
        elif scale == 'logx':
            y2 = splrep(x, np.log(sy), s=0)
            interpolatedData = np.exp(splev(data,y2))
            if plot:
                plt.semilogx(sy, x)
                plt.semilogx(interpolatedData, data, 'bD')
        elif scale == 'logy':
            y2 = splrep(np.log(x), sy, s=0)
            interpolatedData = splev(np.log(data),y2)
            if plot:
                plt.semilogy(sy, x)
                plt.semilogy(interpolatedData, data, 'bD')
        else:
            print(' scale not understood = %s'%(scale))

        if verbose:
            for i, avalue in enumerate(interpolatedData):
                print(' data, value = %12.3e %12.3e'%(data[i], avalue))

        self.IntensityRatioInterpolated = {'data':data, 'value':interpolatedData}


    def gofnt(self, wvlRange=None, top=10, verbose=False, plot = True):
        """
        Calculate the 'so-called' G(T) function.

        Given as a function of both temperature and eDensity.

        Only the top( set by 'top') brightest lines are plotted.
        the G(T) function is returned in a dictionary self.Gofnt

        Note:  if the default value for gui is set to False, then it is usually
            necessary to invoke this method twice to get the desired result.


        Keyword Arguments
        -----------------

        wvlRange:  `array-like`
            the wavelength range to be considered, a two element array-type

        top : `int`
            specifies to plot only the top strongest lines, default = 10

        verbose:  `bool`
            if True, additional information is printed to the console

        plot:  `bool`
            if True, the G(T) functionis plotted, default - True

        """

        if not self.Defaults['gui']:
            print(' it may be necessary to invoke this method twice to ')
            print(' get the desired result')


        if hasattr(self, 'Emiss'):
            em = copy.copy(self.Emiss)
        else:
            self.emiss()
            em = copy.copy(self.Emiss)

        if not hasattr(self, 'Abundance'):
            self.Abundance = io.abundanceRead()

        fontsize = 14
        emiss = em["emiss"]
        wvl = em["wvl"]
        pretty1 = em['pretty1']
        pretty2 = em['pretty2']
        lvl1 = em['lvl1']
        lvl2 = em['lvl2']

        # find which lines are in the wavelength range if it is set
        if wvlRange:
            igvl = util.between(wvl,wvlRange)
        else:
            igvl = range(len(wvl))
        nlines = len(igvl)
        if nlines ==0:
            print(' no lines in selected interval')
            return

        # find the top most intense lines
        if top > nlines:
            top = nlines
        maxEmiss = np.zeros(nlines,np.float64)
        for iline in range(nlines):
            maxEmiss[iline] = emiss[igvl[iline]].max()
        for iline in range(nlines):
            if maxEmiss[iline] >= maxEmiss.max():
                maxAll = emiss[igvl[iline]]

        igvlsort = np.take(igvl,np.argsort(maxEmiss))
        topLines = igvlsort[-top:]
        maxWvl = '%5.3f' % wvl[topLines[-1]]

        # need to make sure there are no negative values before plotting
        good = np.where(emiss > 0.)
        emissMin = emiss[good].min()
        bad = np.where(emiss <= 0.)
        emiss[bad] = emissMin
        topLines = topLines[wvl[topLines].argsort()]
        eDensity = self.EDensity
        temperature = self.Temperature
        ntemp = self.Ntemp

        ndens = self.Ndens

        if plot:
            print(' ndens = %5i ntemp = %5i'%(ndens, ntemp))
            ylabel = 'Emissivity relative to '+maxWvl
            title = self.Spectroscopic

        if ndens==1 and ntemp==1:
            print(' only a single temperature and eDensity')
            return
        elif ndens == 1:
            xlabel = 'Temperature (K)'
            ngofnt = temperature.size
            xvalues = temperature
            outTemperature = temperature
            outDensity = eDensity
            desc_str=' at Density = %10.2e' % eDensity[0]
        elif ntemp == 1:
            xvalues = eDensity
            ngofnt = eDensity.size
            outTemperature = temperature
            outDensity = eDensity
            xlabel = r'$\rm{Electron Density (cm}^{-3}\rm{)}$'
            desc_str = ' at Temperature = %10.2e' % temperature[0]
        else:
            outTemperature = temperature
            outDensity = eDensity
            xlabel = 'Temperature (K)'
            xvalues = temperature
            ngofnt = ntemp
            desc_str = ' for variable Density'

        if plot:
            plt.ion()
            plt.figure()
            ax = plt.subplot(111)
            nxvalues = len(xvalues)
            ymax = 1.2
            ymin = 0.01
            for iline in range(top):
                tline = topLines[iline]
                plt.loglog(xvalues,emiss[tline]/maxAll)
                if np.min(emiss[tline]/maxAll) < ymin:
                    ymin = np.min(emiss[tline]/maxAll)
                skip = 2
                start = divmod(iline,nxvalues)[1]
                for ixvalue in range(start,nxvalues,nxvalues//skip):
                    plt.text(xvalues[ixvalue],emiss[tline,ixvalue]/maxAll[ixvalue],str(wvl[tline]))
            plt.xlim(xvalues.min(),xvalues.max())
            plt.ylim(ymin, ymax)
            plt.xlabel(xlabel,fontsize=fontsize)
            plt.ylabel(ylabel,fontsize=fontsize)
            if ndens == ntemp and ntemp > 1:
                plt.text(0.07, 0.5,title, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
                ax2 = plt.twiny()
                xlabelDen = r'Electron Density (cm$^{-3}$)'
                plt.xlabel(xlabelDen, fontsize=fontsize)
                plt.loglog(eDensity,emiss[topLines[top-1]]/maxAll, visible=False)
                ax2.xaxis.tick_top()
            else:
                plt.ylim(ymin, ymax)
                plt.title(title+desc_str,fontsize=fontsize)
            plt.draw()
            time.sleep(0.5)
        wvlChoices = []
        for iline in range(top):
            tline = topLines[iline]
            wvlChoices.append('%12.4f %4i %4i %s - %s'%(wvl[tline], lvl1[tline], lvl2[tline], pretty1[tline], pretty2[tline]))
#        if plot:
        gline = chGui.gui.selectorDialog(wvlChoices,label='Select line(s)')
        gline_idx = gline.selectedIndex
        #
        #
        gAbund = self.Abundance
        #
        if hasattr(self, 'IoneqOne'):
            thisIoneq = self.IoneqOne
        else:
            self.ioneqOne()
            thisIoneq = self.IoneqOne
        if verbose:
            print(' abundance = %10.2e '%(gAbund))
            print(' index  temperature  ion fraction')
            for it,  anioneq in enumerate(thisIoneq):
                print (' %5i %10.2e %10.2e '%(it, outTemperature[it], anioneq))

        gIoneq = self.IoneqOne/eDensity
        # plot the desired ratio
        if plot:
            plt.figure()
        g_line = topLines[gline_idx]#  [0]
        gofnt = np.zeros(ngofnt,np.float64)
#        if plot:
        for aline in g_line:
            gofnt += gAbund*gIoneq*emiss[aline].squeeze()
        self.Gofnt = {'temperature':outTemperature,'eDensity':outDensity,'gofnt':gofnt, 'index':g_line, 'wvl':wvl[g_line]}
        #
        if plot:
            plt.loglog(xvalues,gofnt)
            plt.xlim(xvalues.min(),xvalues.max())
            plt.xlabel(xlabel,fontsize=fontsize)
            plt.ylabel('Gofnt',fontsize=fontsize)
            newTitle = '%s'%(self.Spectroscopic) + '%12.3f %4i - %4i %s - %s'%(wvl[g_line[0]], lvl1[g_line[0]], lvl2[g_line[0]], pretty1[g_line[0]], pretty2[g_line[0]])
            if len(g_line) > 1:
                newTitle += '\n'
            for igl in g_line[1:]:
                newTitle += ' ' + '%12.3f %4i %4i %s - %s'%(wvl[igl], lvl1[igl], lvl2[igl], pretty1[igl], pretty2[igl])
                if igl != g_line[-1]:
                    newTitle += '\n'
            if ndens == ntemp and ntemp > 1:
                plt.text(0.07, 0.5,newTitle, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
                ax2 = plt.twiny()
                plt.xlabel(xlabelDen, fontsize=fontsize)
                plt.loglog(eDensity,gofnt, visible=False)
                ax2.xaxis.tick_top()
            else:
                plt.title(newTitle, fontsize=fontsize)
            self.Gofnt['transition'] = newTitle
            self.Gofnt['xlabel'] = xlabel
            self.Gofnt['ylabel'] = ylabel
            plt.tight_layout()

    def twoPhotonEmiss(self, wvl):
        """
        To calculate the two-photon continuum rate coefficient - only for hydrogen- and helium-like ions
        """
        if self.Defaults['wavelength'] != 'angstrom':
            print(' the two-photon continuum can only be calculated for wavelengths in angstroms')
            return

        xlabel = self.Labels['xlabel']
        ylabel = self.Labels['emissYlabel']

        if np.array_equal(self.Em, np.ones_like(self.Em)):
            ylabel += '($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'

        wvl = np.array(wvl, np.float64)
        nWvl = wvl.size
        if self.Z -self.Ion > 1 or self.Dielectronic:
            # this is not a hydrogen-like or helium-like ion
            self.TwoPhoton = {'emiss':np.zeros(nWvl, 'float4'), 'wvl':wvl}
            return
        else:
            if hasattr(self, 'Population'):
                pop = self.Population['population']
                nTempDens = max(self.Temperature.size, self.EDensity.size)
            else:
                self.populate()
                pop = self.Population['population']
                nTempDens = max(self.Temperature.size, self.EDensity.size)
            if nTempDens > 1:
                emiss = np.zeros((nTempDens, nWvl), np.float64)
                if self.EDensity.size == 1:
                    eDensity = np.repeat(self.EDensity, nTempDens)
                else:
                    eDensity = self.EDensity
            else:
                emiss = np.zeros(nWvl, np.float64)
                eDensity = self.EDensity
            if self.Z == self.Ion:
                # H seq
                l1 = 1-1
                l2 = 2 - 1
                wvl0 = 1.e+8/(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                goodWvl = wvl > wvl0
                y = wvl0/wvl[goodWvl]
                dist = io.twophotonHRead()
                avalue = dist['avalue'][self.Z-1]
                asum = dist['asum'][self.Z-1]
                distr1 = splrep(dist['y0'], dist['psi0'][self.Z-1], s=0)
                distr = avalue*y*splev(y, distr1)/(asum*wvl[goodWvl])
                if self.Defaults['flux'] == 'energy':
                    f = (const.light*const.planck*1.e+8)/wvl[goodWvl]
                else:
                    f = 1.
                if nTempDens == 1:
                    emiss[goodWvl] = f*pop[l2]*distr/self.EDensity
                else:
                    for it in range(nTempDens):
                        emiss[it, goodWvl] = f*pop[it, l2]*distr/self.EDensity[it]
                self.TwoPhotonEmiss = {'wvl':wvl, 'emiss':emiss}
            else:
                # He seq
                l1 = 1-1
                l2 = heseqLvl2[self.Z -1] - 1
                wvl0 = 1.e+8/(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                goodWvl = wvl > wvl0
                y = wvl0/wvl[goodWvl]
                dist = io.twophotonHeRead()
                avalue = dist['avalue'][self.Z-1]
                distr1 = splrep(dist['y0'], dist['psi0'][self.Z-1], s=0)
                distr = avalue*y*splev(y, distr1)/wvl[goodWvl]
                if self.Defaults['flux'] == 'energy':
                    f = (const.light*const.planck*1.e+8)/wvl[goodWvl]
                else:
                    f = 1.
                for it in range(nTempDens):
                    emiss[it, goodWvl] = f*pop[it, l2]*distr/self.EDensity[it]
                self.TwoPhotonEmiss = {'wvl':wvl, 'emiss':emiss, 'xlabel':xlabel,
                    'ylabel':ylabel}

    def twoPhoton(self, wvl, verbose=False):
        '''
        to calculate the two-photon continuum - only for hydrogen- and helium-like ions
        includes the elemental abundance and the ionization equilibrium
        includes the emission measure if specified
        discussed in [105]_

        References
        ----------

        .. [105] Young et al., 2003, ApJSS, `144, 135  <http://adsabs.harvard.edu/abs/2003ApJS..144..135Y>`_

        '''
        if self.Defaults['wavelength'] != 'angstrom':
            print(' the two-photon continuum can only be calculated for wavelengths in angstroms')
            return

        xlabel = self.Labels['xlabel']
        ylabel = self.Labels['spectrumYlabel']

        if np.array_equal(self.Em, np.ones_like(self.Em)):
            ylabel += '($\int\,$ N$_e\,$N$_H\,$d${\it l}$)$^{-1}$'

        wvl = np.array(wvl, np.float64)
        #
        #
        nWvl = wvl.size
        if self.Z - self.Ion > 1 or self.Dielectronic:
            # this is not a hydrogen-like or helium-like ion
            if verbose:
                print(' not doing 2 photon for %s'%(self.IonStr))
            self.TwoPhoton = {'emiss':np.zeros(nWvl, np.float64), 'wvl':wvl}
            return
        else:
            if hasattr(self, 'Abundance'):
                ab = self.Abundance
            else:
                self.Abundance = io.abundanceRead()
                ab = self.Abundance
            if hasattr(self, 'IoneqOne'):
                thisIoneq = self.IoneqOne
            else:
                self.ioneqOne()
                thisIoneq = self.IoneqOne
            if hasattr(self, 'Population'):
                pop = self.Population['population']
                nTempDens = max(self.Temperature.size, self.EDensity.size)
            else:
                self.populate()
                pop = self.Population['population']
                nTempDens = max(self.Temperature.size, self.EDensity.size)
#            if nTempDens > 1:
            rate = np.zeros((nTempDens, nWvl), np.float64)
            if self.EDensity.size == 1:
                eDensity = np.repeat(self.EDensity, nTempDens)
            else:
                eDensity = self.EDensity
            if self.Temperature.size == 1:
                temperature = np.repeat(self.Temperature, nTempDens)
                thisIoneq = np.repeat(thisIoneq, nTempDens)
            else:
                temperature = self.Temperature
#            else:
#                rate = np.zeros(nWvl, np.float64)
#                eDensity = self.EDensity
#                temperature = self.Temperature
            if self.Z == self.Ion:
                # H seq
                l1 = 1-1
                l2 = 2 - 1
                wvl0 = 1.e+8/(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                goodWvl = wvl > wvl0
                if goodWvl.sum() > 0:
                    y = wvl0/wvl[goodWvl]
                    dist = io.twophotonHRead()
                    avalue = dist['avalue'][self.Z-1]
                    asum = dist['asum'][self.Z-1]
                    distr1 = splrep(dist['y0'], dist['psi0'][self.Z-1], s=0)
                    #  distr is normalized to 2.0
                    distr = avalue*y*splev(y, distr1)/(asum*wvl[goodWvl])
                    if self.Defaults['flux'] == 'energy':
                        f = (const.light*const.planck*1.e+8)/(4.*const.pi*wvl[goodWvl])
                    else:
                        f = 1./(4.*const.pi)
                    for it in range(nTempDens):
                        rate[it, goodWvl] = f*pop[it, l2]*distr*ab*thisIoneq[it]*self.Em[it]/eDensity[it]
                self.TwoPhoton = {'wvl':wvl, 'intensity':rate.squeeze(), 'em':self.Em, 'xlabel':xlabel,
                'ylabel':ylabel}

            else:
                # He seq
                l1 = 1-1
                l2 = heseqLvl2[self.Z -1] - 1
                if verbose:
                    print('two photon l2:  %5i %s  for Z = %5i'%(l2, self.Elvlc['pretty'][l2],  self.Z))
                wvl0 = 1.e+8/(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                goodWvl = wvl > wvl0
                if goodWvl.sum() > 0:
                    y = wvl0/wvl[goodWvl]
                    dist = io.twophotonHeRead()
                    avalue = dist['avalue'][self.Z-1]
                    distr1 = splrep(dist['y0'], dist['psi0'][self.Z-1], s=0)
                    distr = avalue*y*splev(y, distr1)/wvl[goodWvl]
                    if self.Defaults['flux'] == 'energy':
                        f = (const.light*const.planck*1.e+8)/(4.*const.pi*wvl[goodWvl])
                    else:
                        f = 1./(4.*const.pi)
                    for it in range(nTempDens):
                        rate[it, goodWvl] = f*pop[it, l2]*distr*ab*thisIoneq[it]*self.Em[it]/eDensity[it]
                self.TwoPhoton = {'wvl':wvl, 'intensity':rate.squeeze(), 'em':self.Em, 'xlabel':xlabel,
                'ylabel':ylabel}

    def twoPhotonLoss(self):
        '''
        to calculate the two-photon energy loss rate - only for hydrogen- and helium-like ions
        includes the elemental abundance and the ionization equilibrium
        does not include the emission measure
        '''
        if self.Defaults['wavelength'] != 'angstrom':
            print(' the two-photon continuum can only be calculated for wavelengths in angstroms')
            return

        if self.Z -self.Ion > 1 or self.Dielectronic:
            # this is not a hydrogen-like or helium-like ion
            nTempDens = max(self.Temperature.size, self.EDensity.size)
#            if nTempDens > 1:
            rate = np.zeros((nTempDens), np.float64)
            self.TwoPhotonLoss = {'rate':rate}
        else:
            if hasattr(self, 'Abundance'):
                ab = self.Abundance
            else:
                self.Abundance = io.abundanceRead()
                ab = self.Abundance
            if hasattr(self, 'IoneqOne'):
                thisIoneq = self.IoneqOne
            else:
                self.ioneqOne()
                thisIoneq = self.IoneqOne
            if hasattr(self, 'Population'):
                pop = self.Population['population']
                nTempDens = max(self.Temperature.size, self.EDensity.size)
            else:
                self.populate()
                pop = self.Population['population']
                nTempDens = max(self.Temperature.size, self.EDensity.size)
            rate = np.zeros((nTempDens), np.float64)
            eDensity = self.EDensity
#            else:
            eDensity = self.EDensity
            if self.Z == self.Ion:
                # H seq
                l1 = 1 - 1
                l2 = 2 - 1
                wvl0 = 1.e+8/(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                dist = io.twophotonHRead()
                avalue = dist['avalue'][self.Z-1]
                f = (avalue*const.light*const.planck*1.e+8)/wvl0
                for it in range(nTempDens):
                    rate[it] = f*pop[it, l2]*ab*thisIoneq[it]/eDensity[it]
                self.TwoPhotonLoss = {'temperature':self.Temperature,'eDensity':self.EDensity,'rate':rate}
            else:
                # He seq
                l1 = 1 - 1
                l2 = heseqLvl2[self.Z -1] -1
                wvl0 = 1.e+8/(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                dist = io.twophotonHeRead()
                avalue = dist['avalue'][self.Z-1]
                f = (avalue*const.light*const.planck*1.e+8)/wvl0
                for it in range(nTempDens):
                    rate[it] = f*pop[it, l2]*ab*thisIoneq[it]/eDensity[it]
                self.TwoPhotonLoss = {'temperature':self.Temperature,'eDensity':self.EDensity,'rate':rate}
