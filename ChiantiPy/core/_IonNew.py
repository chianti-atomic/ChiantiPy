'''
IonNew
for calculating level populations including dielectronic recombination
using the autoionization rates in the .auto files
'''
import os
import numpy as np
import ChiantiPy.tools.data as chdata
import ChiantiPy.tools.util as util
import ChiantiPy.tools.io as io
    #
    # -------------------------------------------------------------------------
    #
def setupNew(self, dir=0, verbose=0):
    '''
    if ion is initiated with setup=0, this allows the setup to be done at a later point
    perhaps, more importantly,  by setting dir to a directory containing the necessary files
    for a ChiantiPy ion, it allows one to setup an ion with files not in the current
    Chianti directory
    this is a development version for integrating auoionizing levels etc into the main ion
    '''
    #
    # read in all data if in masterlist
    #  if not, there should still be ionization and recombination rates
    #
    MasterList = chdata.MasterList
    #
    if dir:
        fileName = os.path.join(dir, self.IonStr)
    else:
        fileName = util.ion2filename(self.IonStr)
    if self.IonStr in MasterList:
        if dir:
            self.Elvlc = io.elvlcRead('', filename=fileName+'.elvlc',  verbose=verbose)
            self.Wgfa = io.wgfaRead('', filename=fileName+'.wgfa', elvlcname=fileName+'.elvlc', total=1)
            self.Nwgfa=len(self.Wgfa['lvl1'])
            nlvlWgfa = max(self.Wgfa['lvl2'])
            nlvlList =[nlvlWgfa]
#                splupsfile = fileName + '.splups'
            scupsfile = fileName + '.scups'
            if os.path.isfile(scupsfile):
                # happens the case of fe_3 and prob. a few others
                self.Scups = io.scupsRead('', filename=fileName+'.scups')
                self.Nscups=len(self.Scups['lvl1'])
                nlvlScups = max(self.Scups['lvl2'])
                nlvlList.append(nlvlScups)
                self.Nsplups = 0
#                    nlvlSplups = 0
#                else:
#                    if os.path.isfile(splupsfile):
#                        self.Nscups = 0
#                        nlvlScups = 0
#                        # happens the case of fe_3 and prob. a few others
#                        self.Splups = io.splupsRead('', filename=fileName+'.splups')
#                        self.Nsplups=len(self.Splups['lvl1'])
#                        nlvlSplups = max(self.Splups['lvl2'])
#                        nlvlList.append(nlvlSplups)
            else:
                self.Nscups = 0
                nlvlScups = 0
                print('do not have a scups file for %s'%(self.IonStr))
        else:
            self.Elvlc = io.elvlcRead(self.IonStr,  verbose=verbose)
            self.Wgfa = io.wgfaRead(self.IonStr, total=1)
            self.Nwgfa=len(self.Wgfa['lvl1'])
            nlvlWgfa = max(self.Wgfa['lvl2'])
            nlvlList =[nlvlWgfa]
#                splupsfile = fileName + '.splups'
            scupsfile = fileName + '.scups'
            if os.path.isfile(scupsfile):
                # happens the case of fe_3 and prob. a few others
                self.Scups = io.scupsRead(self.IonStr)
                self.Nsplups=len(self.Scups['lvl1'])
                nlvlScups = max(self.Scups['lvl2'])
                nlvlList.append(nlvlScups)
                self.Nsplups = 0
#                    nlvlSplups = 0
#                else:
#                    if os.path.isfile(splupsfile):
#                        self.Nscups = 0
#                        nlvlScups = 0
#                        # happens the case of fe_3 and prob. a few others
#                        self.Splups = io.splupsRead(self.IonStr)
#                        self.Nsplups=len(self.Splups['lvl1'])
#                        nlvlSplups = max(self.Splups['lvl2'])
#                        nlvlList.append(nlvlSplups)
            else:
                print('do not have either a scups  file for %s'(self.IonStr))
                self.Nscups = 0
                nlvlScups = 0

##                self.Nlvls = nlvlElvlc
        #
        file = fileName +'.cilvl'
        if os.path.isfile(file):
            self.Cilvl = io.cireclvlRead('',filename = fileName, cilvl=1)
            self.Ncilvl=len(self.Cilvl['lvl1'])
            nlvlCilvl = max(self.Cilvl['lvl2'])
            nlvlList.append(nlvlCilvl)
        else:
            self.Ncilvl = 0
        #
        #  not longer using the reclvl files
        #  using the new rrlvl - radiaitive recombination rates only
        #  dielectronic rates derived from the autoionization values in the .auto file
        reclvlfile = fileName +'.reclvl'
        if os.path.isfile(reclvlfile):
            self.Reclvl = io.cireclvlRead('',filename=fileName, reclvl=1)
            self.Nreclvl = len(self.Reclvl['lvl1'])
            nlvlReclvl = max(self.Reclvl['lvl2'])
            nlvlList.append(nlvlReclvl)
        else:
            self.Nreclvl = 0
        # in setupNew, the dielsplups files are disregarded
        #  .dielsplups file may not exist
#            dielsplupsfile = fileName +'.dielsplups'
#            if os.path.isfile(dielsplupsfile):
#                self.DielSplups = io.splupsRead('', filename=dielsplupsfile, filetype='splups')
#                self.Ndielsplups=len(self.DielSplups["lvl1"])
#                nlvlDielSplups = max(self.DielSplups['lvl2'])
#                nlvlList.append(nlvlDielSplups)
#            else:
#                self.Ndielsplups = 0
        #
        #  psplups file may not exist
        psplupsfile = fileName +'.psplups'
        if os.path.isfile(psplupsfile):
            self.Psplups = io.splupsRead('', filename=psplupsfile,  filetype='psplups')
            self.Npsplups=len(self.Psplups["lvl1"])
        else:
            self.Npsplups = 0
        #
        drparamsFile = fileName +'.drparams'
        if os.path.isfile(drparamsFile):
            self.DrParams = io.drRead(self.IonStr)
        #
        rrparamsFile = fileName +'.rrparams'
        if os.path.isfile(rrparamsFile):
            self.RrParams = io.rrRead(self.IonStr)
        #
        # get autoionizing A-values
        autoFile = fileName + '.auto'
        if os.path.isfile(autoFile):
            self.Auto = io.wgfaRead('', filename=autoFile)

        #  not needed for ion, only phion
#                photoxfile = util.ion2filename(self.IonStr)+'.photox'
#                if os.path.isfile(photoxfile):
#                    self.Photox = util.photoxRead(self.IonStr)
        #
        # need to determine the number of levels that can be populated
        nlvlElvlc = len(self.Elvlc['lvl'])
#                print ' nlvlElvlc = ', nlvlElvlc
#                print ' other nlvls = ',  nlvlList
#                nlvlWgfa = max(self.Wgfa['lvl2'])
        #  elvlc file can have more levels than the rate level files
        self.Nlvls = min([nlvlElvlc, max(nlvlList)])
    #
    # -------------------------------------------------------------------------
    #
def populateNew(self, popCorrect=1, verbose=0, **kwargs):
    """
    Calculate level populations for specified ion.  This is a new version that will enable the calculation
    of dielectronic satellite lines without resorting to the dielectronic ions, such as c_5d
    possible keyword arguments include temperature, eDensity, pDensity, radTemperature and rStar
    this is a developmental method for using the autoionizing A-values to determine level resolved
    dielectronic recombination rates
    """
    #
    #
    for one in kwargs.keys():
        if one not in chdata.keywordArgs:
            print(' keyword is not understood - %s'%(one))
    #
    nlvls = self.Nlvls
    nwgfa = self.Nwgfa
    nscups = self.Nscups
    npsplups = self.Npsplups
    #
    if kwargs.has_key('temperature'):
        self.Temperature = np.asarray(kwargs['temperature'])
        temperature = self.Temperature
    elif hasattr(self, 'Temperature'):
        temperature=self.Temperature
    else:
        print(' no temperature values have been set')
        return {'errorMessage':' no temperature values have been set'}
    #
    if kwargs.has_key('eDensity'):
        self.EDensity = np.asarray(kwargs['eDensity'])
        eDensity = self.EDensity
    elif hasattr(self, 'EDensity'):
        eDensity = self.EDensity
    else:
        print(' no eDensity values have been set')
        return {'errorMessage':' no eDensity values have been set'}
    #
    if kwargs.has_key('pDensity'):
        if kwargs['pDensity'] == 'default':
            self.p2eRatio()
            protonDensity = self.ProtonDensityRatio*self.EDensity
        else:
            try:
                self.PDensity = np.asarray(kwargs['pDensity'])
            except:
                print(' could not interpret value for keyword pDensity')
                print(' should be either "default" or a number or array')
                return
    else:
        if hasattr(self, 'PDensity'):
            protonDensity = self.PDensity
        else:
            self.p2eRatio()
            self.PDensity = self.ProtonDensityRatio*self.EDensity
            protonDensity = self.PDensity
            print(' proton density not specified, set to "default"')
    #
    if 'radTemperature' in kwargs.keys() and 'rStar' in kwargs.keys():
        self.RadTemperature = np.asarray(kwargs['radTemperature'])
        radTemperature = np.array(self.RadTemperature)
        self.RStar = np.asarray(kwargs['rStar'])
        rStar = np.asarray(self.RStar)
    elif hasattr(self, 'RadTemperature') and hasattr(self, 'RStar'):
        radTemperature = self.RadTemperature
        rStar = self.RStar
    #
    #
    #
    if self.Ncilvl:
        ci = 1
        cilvl = self.Cilvl
#            if hasattr(self, 'CilvlRate'):
#                cilvlRate = self.CilvlRate
#            else:
#                self.cireclvlDescale('cilvl')
#                cilvlRate = self.CilvlRate
        self.recombRate()
        #
        lowers = util.zion2name(self.Z, self.Ion-1)
        # get the lower ionization stage
        self.Lower = ion(lowers, temperature=self.Temperature, eDensity = self.EDensity)
        self.Lower.ionizRate()
        # need to get multiplicity of lower ionization stage
        lowMult = self.Lower.Elvlc['mult']
    else:
        ci = 0
    #  evetually will be looking for just an recLvl attribute
    #
    #  if the higher ion does not exist in the database, rec=0
    highers = util.zion2name(self.Z, self.Ion+1)
    if self.Nreclvl:
        reclvl = self.Reclvl
        if hasattr(self, 'ReclvlRate'):
            reclvlRate = self.ReclvlRate
        else:
            self.cireclvlDescale('reclvl')
            reclvlRate = self.ReclvlRate
    if hasattr(self, 'Auto'):
        self.drRateLvl()
    if self.Nreclvl or hasattr(self, 'Auto'):
        rec=1
        # get ionization rate of this ion
        self.ionizRate()
            #  get the higher ionization stage
    else:
        rec = 0
    #  if the higher ion does not exist in the database, rec=0
    highers = util.zion2name(self.Z, self.Ion+1)
#        if highers in chdata.MasterList:
    self.Higher = ion(highers, temperature=self.Temperature, eDensity=self.EDensity)
    self.Higher.recombRate()
#        else:
#        #
    rad=np.zeros((nlvls+ci+rec,nlvls+ci+rec),"float64")  #  the populating matrix for radiative transitions
    #
    #
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
                if self.Elvlc['ecm'][l2] >= 0.:
                    ecm2 = self.Elvlc['ecm'][l2]
                else:
                    ecm2 = self.Elvlc['ecmth'][l2]
                #
                if self.Elvlc['ecm'][l1] >= 0.:
                    ecm1 = self.Elvlc['ecm'][l1]
                else:
                    ecm1 = self.Elvlc['ecmth'][l1]
                de = const.invCm2Erg*(ecm2 - ecm1)
                dekt = de/(const.boltzmann*self.RadTemperature)
                # photoexcitation
                phexFactor = dilute*(float(self.Elvlc['mult'][l2])/float(self.Elvlc['mult'][l1]))/(np.exp(dekt) -1.)
                rad[l2+ci,l1+ci] += self.Wgfa["avalue"][iwgfa]*phexFactor
                rad[l1+ci,l1+ci] -= self.Wgfa["avalue"][iwgfa]*phexFactor
                # stimulated emission
                stemFactor = dilute/(np.exp(-dekt) -1.)
                rad[l1+ci,l2+ci] += self.Wgfa["avalue"][iwgfa]*stemFactor
                rad[l2+ci,l2+ci] -= self.Wgfa["avalue"][iwgfa]*stemFactor
    if hasattr(self, 'Auto'):
        # as of 7/2012, we only consider the ground level in the next higher ionization stage
        # hence, requiring lvl1 = 1 or, l1 = 0
        for iauto, avalue in enumerate(self.Auto['avalue']):
            l1 = self.Auto["lvl1"][iauto] -1
            l2 = self.Auto["lvl2"][iauto] -1
            # for now only consider a single level for upper/higher ion
            if l1 == 0 and rec:
                rad[l1 + ci + nlvls, l2 + ci] += avalue
                rad[l2 + ci,  l2 + ci] -= avalue
            # if the higher ion is not in the database, decay to the ground level
            elif l1 == 0:
                rad[l1 + ci, l2 + ci] += avalue
                rad[l2 + ci,  l2 + ci] -= avalue

    #
    #
    if self.Nscups:
#            print(' Nscups = %10i'%(self.Nscups))
        self.upsilonDescale()
        ups = self.Upsilon['upsilon']
        exRate = self.Upsilon['exRate']
        dexRate = self.Upsilon['dexRate']
    #
    if self.Npsplups:
        self.upsilonDescaleSplups(prot=1)
#            pups = self.PUpsilon['upsilon']
        pexRate = self.PUpsilon['exRate']
        pdexRate = self.PUpsilon['dexRate']
    #
    temp=temperature
    ntemp=temp.size
    #
    cc=const.collision*self.EDensity
    ndens=cc.size
    if self.Npsplups:
        cp=const.collision*protonDensity
    if ntemp > 1 and ndens >1 and ntemp != ndens:
        print(' unless temperature or eDensity are single values')
        print(' the number of temperatures values must match the ')
        print(' the number of eDensity values')
        return
    #
    # get corrections for recombination and excitation
    #
    nscups = self.Nscups
    #
    #  first, for ntemp=ndens=1
    if ndens==1 and ntemp==1:
        popmat=np.copy(rad)
        for iscups in range(0,nscups):
            l1=self.Scups["lvl1"][iscups]-1
            l2=self.Scups["lvl2"][iscups]-1
            #
            popmat[l1+ci,l2+ci] += self.EDensity*dexRate[iscups]
            popmat[l2+ci,l1+ci] += self.EDensity*exRate[iscups]
            popmat[l1+ci,l1+ci] -= self.EDensity*exRate[iscups]
            popmat[l2+ci,l2+ci] -= self.EDensity*dexRate[iscups]
            #
        for isplups in range(0,npsplups):
            l1=self.Psplups["lvl1"][isplups]-1
            l2=self.Psplups["lvl2"][isplups]-1
             #
            popmat[l1+ci,l2+ci] += self.PDensity*pdexRate[isplups]
            popmat[l2+ci,l1+ci] += self.PDensity*pexRate[isplups]
            popmat[l1+ci,l1+ci] -= self.PDensity*pexRate[isplups]
            popmat[l2+ci,l2+ci] -= self.PDensity*pdexRate[isplups]
       # now include ionization rate from
        if ci:
            #
            # the ciRate can be computed for all temperatures
            #
            ciTot = 0.
            for itrans in range(len(cilvl['lvl1'])):
                lvl1 = cilvl['lvl1'][itrans]-1
                lvl2 = cilvl['lvl2'][itrans]-1
                # this is kind of double booking the ionization rate components
                popmat[lvl2+ci, lvl1] += self.EDensity*self.CilvlRate['rate'][itrans]
                popmat[lvl1, lvl1] -= self.EDensity*self.CilvlRate['rate'][itrans]
                ciTot += self.EDensity*self.CilvlRate['rate'][itrans]
            #
            popmat[1, 0] += (self.EDensity*self.Lower.IonizRate['rate'] - ciTot)
            popmat[0, 0] -= (self.EDensity*self.Lower.IonizRate['rate'] - ciTot)
            popmat[0, 1] += self.EDensity*self.RecombRate['rate']
            popmat[1, 1] -= self.EDensity*self.RecombRate['rate']
        if rec:
            # ntemp=ndens=1
            #
            if hasattr(self, 'DrRateLvl'):
#                    branch = np.zeros(self.Ndielsplups, 'float64')
                for idr, rate in enumerate(self.DrRateLvl['rate']):
                    l1 = self.Auto["lvl1"][idr] - 1
                    l2 = self.Auto["lvl2"][idr] - 1
#                        popmat[l2+ci,-1] += self.EDensity*self.DrRateLvl['rate'][idr]
#                        popmat[-1, -1] -= self.EDensity*self.DrRateLvl['rate'][idr]
                    popmat[l2+ci,-1] += self.EDensity*rate
                    popmat[-1, -1] -= self.EDensity*rate
                    #
                dielTot = self.DrRateLvl['totalRate'][0]
            else:
                dielTot = 0.
            #
            #
            #
            for itrans in range(self.Nreclvl):
#                    lvl1 = reclvl['lvl1'][itrans]-1
                lvl2 = reclvl['lvl2'][itrans]-1
                popmat[lvl2+ci, -1] += self.EDensity*reclvlRate['rate'][itrans]
                popmat[-1, -1] -= self.EDensity*reclvlRate['rate'][itrans]
            if self.Nreclvl:
                recTot = reclvlRate['rate'].sum(axis=0)
            else:
                recTot = 0.
            #
            #
            popmat[-1,  ci] += self.EDensity*self.IonizRate['rate']
            popmat[ci, ci] -= self.EDensity*self.IonizRate['rate']
            #
            # next 2 line take care of overbooking
            netRecomb = self.EDensity*(self.Higher.RecombRate['rate'][0]- recTot - dielTot)
            #
            if netRecomb > 0.:
                popmat[ci, -1] += netRecomb
                popmat[-1, -1] -= netRecomb
        #
        # normalize to unity
        norm=np.ones(nlvls+ci+rec,'float64')
        if ci:
            norm[0] = 0.
        if rec:
            norm[nlvls+ci+rec-1] = 0.
        popmata = np.copy(popmat)
        normRow = (nlvls+ci+rec-1)/2
        #popmata[nlvls+ci+rec-1]=norm
        popmata[normRow]=norm
        b=np.zeros(nlvls+ci+rec,'float64')
        #b[nlvls+ci+rec-1]=1.
        b[normRow] = 1.
        try:
            fullpop=np.linalg.solve(popmata,b)
            pop = fullpop[ci:ci+nlvls]
            if rec:
                popHigher = fullpop[-1]
            else:
                popHigher = 0.
        except np.linalg.LinAlgError:
            pop = np.zeros(nlvls, 'float64')
            popHigher = 0.
#                print ' error in matrix inversion, setting populations to zero at T = ', ('%8.2e')%(temperature)
    #
    # ------------- ntemp = 1 ---------------------------------------------------------
    #
    #
    elif ndens == 1:
        pop = np.zeros((ntemp, nlvls),"float64")
        popHigher = np.zeros(ntemp, 'float64')
#            pop=np.zeros((ntemp,ci + nlvls + rec),"float64")
        for itemp in range(ntemp):
            popmat=np.copy(rad)
            for iscups in range(0,nscups):
                l1=self.Scups["lvl1"][iscups]-1
                l2=self.Scups["lvl2"][iscups]-1
                popmat[l1+ci,l2+ci] += self.EDensity*dexRate[iscups, itemp]
                popmat[l2+ci,l1+ci] += self.EDensity*exRate[iscups, itemp]
                popmat[l1+ci,l1+ci] -= self.EDensity*exRate[iscups, itemp]
                popmat[l2+ci,l2+ci] -= self.EDensity*dexRate[iscups, itemp]
            for isplups in range(0,npsplups):
                l1=self.Psplups["lvl1"][isplups]-1
                l2=self.Psplups["lvl2"][isplups]-1
                # for proton excitation, the levels are all below the ionization potential
                 #
                popmat[l1+ci,l2+ci] += self.PDensity[itemp]*pdexRate[isplups, itemp]
                popmat[l2+ci,l1+ci] += self.PDensity[itemp]*pexRate[isplups, itemp]
                popmat[l1+ci,l1+ci] -= self.PDensity[itemp]*pexRate[isplups, itemp]
                popmat[l2+ci,l2+ci] -= self.PDensity[itemp]*pdexRate[isplups, itemp]
            # now include ionization rate from
            if ci:
                #
                # the ciRate can be computed for all temperatures
                #
                ciTot = 0.
                for itrans in range(len(cilvl['lvl1'])):
                    lvl1 = cilvl['lvl1'][itrans]-1
                    lvl2 = cilvl['lvl2'][itrans]-1
                    #mult = lowMult[lvl1-1]
                    popmat[lvl2+ci, lvl1] += self.EDensity*self.CilvlRate['rate'][itrans, itemp]
                    popmat[lvl1, lvl1] -= self.EDensity*self.CilvlRate['rate'][itrans, itemp]
                    ciTot += self.EDensity*self.CilvlRate['rate'][itrans, itemp]
                #
                popmat[1, 0] += (self.EDensity*self.Lower.IonizRate['rate'][itemp] - ciTot)
                popmat[0, 0] -= (self.EDensity*self.Lower.IonizRate['rate'][itemp] - ciTot)
                popmat[0, 1] += self.EDensity*self.RecombRate['rate'][itemp]
                popmat[1, 1] -= self.EDensity*self.RecombRate['rate'][itemp]
            if rec:
            #
            # ndens=1
                if hasattr(self, 'DrRateLvl'):
                    for idr, rate in enumerate(self.DrRateLvl['rate']):
                        l1 = self.Auto["lvl1"][idr] - 1
                        l2 = self.Auto["lvl2"][idr] - 1
                        popmat[l2 + ci, -1] += self.EDensity*rate[itemp]
                        popmat[-1, -1] -= self.EDensity*rate[itemp]
                        #
                    # DrRateLvl['totalRate'] includes the branching ratio
                    dielTot = self.DrRateLvl['totalRate'][itemp]
                else:
                    dielTot = 0.
                #
                #
                for itrans in range(self.Nreclvl):
                    lvl1 = reclvl['lvl1'][itrans]-1
                    lvl2 = reclvl['lvl2'][itrans]-1
                    popmat[lvl2+ci, -1] += self.EDensity*self.ReclvlRate['rate'][itrans, itemp]
                    popmat[-1, -1] -= self.EDensity*self.ReclvlRate['rate'][itrans, itemp]
                #
                if self.Nreclvl:
                    recTot = self.ReclvlRate['rate'][:, itemp].sum()
                else:
                    recTot = 0.
                #
                #
                popmat[-1, ci] += self.EDensity*self.IonizRate['rate'][itemp]
                popmat[ci, ci] -= self.EDensity*self.IonizRate['rate'][itemp]
                #
                netRecomb = self.EDensity*(self.Higher.RecombRate['rate'][itemp]- recTot - dielTot)
                #
                if netRecomb > 0.:
                    popmat[ci, -1] += netRecomb
                    popmat[-1, -1] -= netRecomb
                #
            # normalize to unity
            # ndens = 1
            norm=np.ones(nlvls+ci+rec,'float64')
            if ci:
                norm[0] = 0.
            if rec:
                norm[-1] = 0.
            popmata = np.copy(popmat)
            normRow = (nlvls+ci+rec-1)/2
            #popmata[nlvls+ci+rec-1]=norm
            popmata[normRow]=norm
            b=np.zeros(nlvls+ci+rec,'float64')
            #b[nlvls+ci+rec-1]=1.
            b[normRow] = 1.
            #popmata[nlvls+ci+rec-1]=norm
            #b=np.zeros(nlvls+ci+rec,'float64')
            #b[nlvls+ci+rec-1] = 1.
#                b[-1] = 1.
            try:
                thispop=np.linalg.solve(popmata,b)
                pop[itemp] = thispop[ci:ci+nlvls]
                popHigher[itemp] = thispop[-1]
            except np.linalg.LinAlgError:
                pop[itemp] = np.zeros(nlvls, 'float64')
                popHigher[itemp] = 0.
#                    print ' error in matrix inversion, setting populations to zero at T = ', ('%8.2e')%(temperature[itemp])
        #
    elif ntemp == 1:
        pop=np.zeros((ndens,nlvls),"float64")
        popHigher = np.zeros(ndens, 'float64')
        for idens in range(0,ndens):
            popmat=np.copy(rad)
            for iscups in range(0,nscups):
                l1=self.Scups["lvl1"][iscups]-1
                l2=self.Scups["lvl2"][iscups]-1
            #
                popmat[l1+ci,l2+ci] += self.EDensity[idens]*dexRate[iscups]
                popmat[l2+ci,l1+ci] += self.EDensity[idens]*exRate[iscups]
                popmat[l1+ci,l1+ci] -= self.EDensity[idens]*exRate[iscups]
                popmat[l2+ci,l2+ci] -= self.EDensity[idens]*dexRate[iscups]
            #
            for isplups in range(0,npsplups):
                l1=self.Psplups["lvl1"][isplups]-1
                l2=self.Psplups["lvl2"][isplups]-1
             #
                popmat[l1+ci,l2+ci] += self.PDensity[idens]*pdexRate[isplups]
                popmat[l2+ci,l1+ci] += self.PDensity[idens]*pexRate[isplups]
                popmat[l1+ci,l1+ci] -= self.PDensity[idens]*pexRate[isplups]
                popmat[l2+ci,l2+ci] -= self.PDensity[idens]*pdexRate[isplups]
            # now include ionization rate from
            if ci:
                #
                #
                ciTot = 0.
                for itrans in range(len(cilvl['lvl1'])):
                    lvl1 = cilvl['lvl1'][itrans] -1
                    lvl2 = cilvl['lvl2'][itrans] -1
                    popmat[lvl2+ci, lvl1] += self.EDensity[idens]*self.CilvlRate['rate'][itrans]
                    popmat[lvl1, lvl1] -= self.EDensity[idens]*self.CilvlRate['rate'][itrans]
                    ciTot += self.EDensity[idens]*self.CilvlRate['rate'][itrans]
                popmat[1, 0] += (self.EDensity[idens]*self.Lower.IonizRate['rate'] -ciTot)
                popmat[0, 0] -= (self.EDensity[idens]*self.Lower.IonizRate['rate'] -ciTot)
                popmat[0, 1] += self.EDensity[idens]*self.RecombRate['rate']
                popmat[1, 1] -= self.EDensity[idens]*self.RecombRate['rate']
            if rec:
                #ntemp = 1
                if hasattr(self, 'DrRateLvl'):
                    for idr, rate in enumerate(self.DrRateLvl['rate']):
                        l1 = self.Auto["lvl1"][idr] - 1
                        l2 = self.Auto["lvl2"][idr] - 1
#                            popmat[l2+ci,l1+ci+nlvls] += self.EDensity[idens]*self.DrRateLvl['rate'][idr]
#                            popmat[l1+ci,l1+ci] -= self.EDensity[idens]*self.DrRateLvl['rate'][idr]
                        popmat[l2+ci,-1] += self.EDensity[idens]*rate
                        popmat[-1, -1] -= self.EDensity[idens]*rate
                        #
                    dielTot = self.DrRateLvl['totalRate'][0]
                else:
                    dielTot = 0.
                if self.Nreclvl:
                    recTot = self.ReclvlRate['rate'].sum()
                else:
                    recTot = 0.
                #
                popmat[-1,  ci] += self.EDensity[idens]*self.IonizRate['rate']
                popmat[ci, ci] -= self.EDensity[idens]*self.IonizRate['rate']
                #
                netRecomb = self.EDensity[idens]*(self.Higher.RecombRate['rate'] - recTot - dielTot)
                if netRecomb > 0.:
                    popmat[ci, -1] += netRecomb
                    popmat[-1, -1] -= netRecomb
                #
#                    for itrans in range(len(rrlvl['lvl1'])):
                for itrans in range(self.Nreclvl):
                    lvl1 = reclvl['lvl1'][itrans]-1
                    lvl2 = reclvl['lvl2'][itrans]-1
                    popmat[lvl2+ci, -1] += self.EDensity[idens]*self.ReclvlRate['rate'][itrans]
                    popmat[-1, -1] -= self.EDensity[idens]*self.ReclvlRate['rate'][itrans]
            #
            # normalize to unity
            norm=np.ones(nlvls+ci+rec,'float64')
            if ci:
                norm[0] = 0.
            if rec:
                norm[-1] = 0.
            popmata = np.copy(popmat)
            popmata[nlvls+ci+rec-1] = norm
#                popmat[nlvls+ci+rec-1] = norm
            b = np.zeros(nlvls+ci+rec,'float64')
            b[nlvls+ci+rec-1] = 1.
            try:
                thispop=np.linalg.solve(popmata,b)
                pop[idens] = thispop[ci:ci+nlvls]
                popHigher[idens] = thispop[-1]
            except np.linalg.LinAlgError:
                pop[idens] = np.zeros(nlvls, 'float64')
                popHigher[idens] = 0.
#                    print ' error in matrix inversion, setting populations to zero at eDensity = ', ('%8.2e')%(eDensity[idens])
            #
    elif ntemp>1  and ntemp==ndens:
        pop=np.zeros((ntemp,nlvls),"float64")
        popHigher = np.zeros(ntemp, 'float64')
        for itemp in range(0,ntemp):
            temp=self.Temperature[itemp]
            popmat=np.copy(rad)
            for iscups in range(0,nscups):
                l1=self.Scups["lvl1"][iscups]-1
                l2=self.Scups["lvl2"][iscups]-1
                #
                popmat[l1+ci,l2+ci] += self.EDensity[itemp]*dexRate[iscups, itemp]
                popmat[l2+ci,l1+ci] += self.EDensity[itemp]*exRate[iscups, itemp]
                popmat[l1+ci,l1+ci] -= self.EDensity[itemp]*exRate[iscups, itemp]
                popmat[l2+ci,l2+ci] -= self.EDensity[itemp]*dexRate[iscups, itemp]
            # proton rates
            for isplups in range(0,npsplups):
                l1=self.Psplups["lvl1"][isplups]-1
                l2=self.Psplups["lvl2"][isplups]-1
                # for proton excitation, the levels are all below the ionization potential
                 #
                popmat[l1+ci,l2+ci] += self.PDensity[itemp]*pdexRate[isplups, itemp]
                popmat[l2+ci,l1+ci] += self.PDensity[itemp]*pexRate[isplups, itemp]
                popmat[l1+ci,l1+ci] -= self.PDensity[itemp]*pexRate[isplups, itemp]
                popmat[l2+ci,l2+ci] -= self.PDensity[itemp]*pdexRate[isplups, itemp]
            # now include ionization rate from
            if ci:
                #
                # the ciRate can be computed for all temperatures
                #
                ciTot = 0.
                for itrans in range(len(cilvl['lvl1'])):
                    lvl1 = cilvl['lvl1'][itrans] -1
                    lvl2 = cilvl['lvl2'][itrans] -1
                    popmat[lvl2+ci, lvl1] += self.EDensity[itemp]*self.CilvlRate['rate'][itrans, itemp]
                    popmat[lvl1, lvl1] -= self.EDensity[itemp]*self.CilvlRate['rate'][itrans, itemp]
                    ciTot += self.EDensity[itemp]*self.CilvlRAte['rate'][itrans, itemp]
                popmat[1, 0] += (self.EDensity[itemp]*self.Lower.IonizRate['rate'][itemp] - ciTot)
                popmat[0, 0] -= (self.EDensity[itemp]*self.Lower.IonizRate['rate'][itemp] - ciTot)
                popmat[0, 1] += self.EDensity[itemp]*self.RecombRate['rate'][itemp]
                popmat[1, 1] -= self.EDensity[itemp]*self.RecombRate['rate'][itemp]
            if rec:
            #
                if hasattr(self, 'DrRateLvl'):
#                    branch = np.zeros(self.Ndielsplups, 'float64')
                    for idr, rate in enumerate(self.DrRateLvl['rate']):
                        l1 = self.Auto["lvl1"][idr] - 1
                        l2 = self.Auto["lvl2"][idr] - 1
                        popmat[l2+ci,l1+ci+nlvls] += self.EDensity[itemp]*rate[itemp]
                        popmat[-1, -1] -= self.EDensity[itemp]*rate[itemp]
                        #
                    dielTot = self.DrRateLvl['totalRate'][itemp]
                else:
                    dielTot = 0.
            #
                if self.Nreclvl:
                    recTot = self.RrlvlRate['rate'][:, itemp].sum()
                else:
                    recTot = 0.
            #
                #
                popmat[-1,  ci] += self.EDensity[itemp]*self.IonizRate['rate'][itemp]
                popmat[ci, ci] -= self.EDensity[itemp]*self.IonizRate['rate'][itemp]
                #
                netRecomb = self.EDensity[itemp]*(self.Higher.RecombRate['rate'][itemp] - recTot - dielTot)
                if netRecomb > 0.:
                    popmat[ci, -1] += netRecomb
                    popmat[-1, -1] -= netRecomb
                #
                for itrans in range(self.Nreclvl):
                    lvl1 = reclvl['lvl1'][itrans]-1
                    lvl2 = reclvl['lvl2'][itrans]-1
                    popmat[lvl2+ci, -1] += self.EDensity[itemp]*self.RrlvlRate['rate'][itrans, itemp]
                    popmat[-1, -1] -= self.EDensity[itemp]*self.RrlvlRate['rate'][itrans, itemp]
            #
            # normalize to unity
            norm=np.ones(nlvls+ci+rec,'float64')
            if ci:
                norm[0] = 0.
            if rec:
                norm[-1] = 0.
            popmata = np.copy(popmat)
            popmata[nlvls+ci+rec-1] = norm
#                popmat[nlvls+ci+rec-1] = norm
            b=np.zeros(nlvls+ci+rec,'float64')
            b[nlvls+ci+rec-1]=1.
            try:
                thispop=np.linalg.solve(popmata,b)
                pop[itemp] = thispop[ci:ci+nlvls]
                popHigher[itemp] = thispop[-1]
            except np.linalg.LinAlgError:
                pop[itemp] = np.zeros(nlvls, 'float64')
                popHigher[itemp] = 0.
#                    print ' error in matrix inversion, setting populations to zero at T = ', ('%8.2e')%(temperature[itemp])
        #
    pop=np.where(pop >0., pop,0.)
    self.Population={"temperature":temperature, "eDensity":eDensity, "population":pop, "protonDensity":protonDensity, "ci":ci, "rec":rec, 'popmat':popmata, 'b':b, 'rad':rad}
    if rec:
        self.Population['popHigher']= popHigher
        self.Population['higher'] = self.Higher
        self.Population['recTot'] = recTot
        self.Population['dielTot'] = dielTot
        self.Population['netRecomb'] = netRecomb
    #
    return
    #
    # -------------------------------------------------------------------------------------
    #
