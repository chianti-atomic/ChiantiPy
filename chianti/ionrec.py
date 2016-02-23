#import types
import numpy as np
import chianti.core as ch
from chianti import util
class ioneq(ch.ion):
    def __init__(self,z, temperature, verbose=False):
        '''calculates the ionization equilibrium for element z at a single temperatures'''
#        self.Defaults=defaults
        self.Z=z
        self.Temperature=temperature
        #
    def getIoneq(self):
        ''' get the ionization equilibrium for this ion at the specified temperature
        getIoneqFromEvolver does the same but is a bit less accurate'''
        z = self.Z
        temperature = self.Temperature
        ionList=[]
        chIons=[]
        for ion in range(1, z+2):
            ions=util.zion2name(z, ion)
            ionList.append(ions)
#            print z, ion, ions
            atom=ch.ion(ions, temperature=temperature, setup=0)
            atom.ionizRate()
            atom.recombRate()
            chIons.append(atom)
        ioneq=np.zeros((z+1), 'float64')
        factor=[]
        for anIon in chIons:
            if type(anIon.IonizRate) != type(None) and type(anIon.RecombRate) != type(None):
                rat=anIon.IonizRate['rate']/anIon.RecombRate['rate']
                factor.append(rat**2 + rat**(-2))
            else:
                factor.append(0.)
        factor[0]=max(factor)
        factor[-1]=max(factor)
        ionmax=factor.index(min(factor))
        ioneq[ionmax]=1.
        #
        for iz in range(ionmax+1, z+1):
            ionrate=chIons[iz-1].IonizRate['rate']
            recrate=chIons[iz].RecombRate['rate']
            ioneq[iz]=ionrate*ioneq[iz-1]/recrate
#            print ' iz, ioneq = ',iz,ioneq[iz]
        #
        for iz in range(ionmax-1, -1, -1):
            ionrate=chIons[iz].IonizRate['rate']
            recrate=chIons[iz+1].RecombRate['rate']
            ioneq[iz]=recrate*ioneq[iz+1]/ionrate
#            print ' iz, ioneq = ',iz,ioneq[iz]
        ionsum=ioneq.sum()
#        print ' ionsum = ', ionsum
        ioneq=ioneq/ionsum
        self.Ioneq=ioneq
        #
        # --------------------------------------
        #
    def getIoneqFromStatic(self):
        '''gets the ionization equilibrium for the specified ion at the specified temperature
        not as accurate as getIoneq since this methods uses a matrix inversion which can result in
        nonphysical values (usually small) that are zero'd out '''
        z = self.Z
        self.getEvolver(self.Temperature)
        b = np.zeros((z+1),'float64')
        b[0] = 1.
        static = self.Static
        ioneq = np.linalg.solve(evolver,b)
        ioneq = np.where(ioneq >= 0., ioneq, 0.)
        self.Ioneq = ioneq
        #
        # --------------------------------------
        #

    def getStatic(self):
        ''' calculate the static matrix for the specified temperature'''
        temperature = self.Temperature
        z = self.Z
        ionList=[]
        chIons=[]
        evolver = np.zeros((z+1,z+1),'float64')
        for ion in range(1, z+2):
            ions=util.zion2name(z, ion)
            ionList.append(ions)
#            print z, ion, ions
            atom=ch.ion(ions, temperature=temperature)
            atom.ionizRate()
            atom.recombRate()
            chIons.append(atom)
        for stage in range(0,z):
            atom = chIons[stage]
#            if type(atom.IonizRate) != types.NoneType:
            evolver[stage,stage] -= atom.IonizRate['rate']
            atom = chIons[stage+1]
            evolver[stage,stage+1] += atom.RecombRate['rate']
#                evolver[stage-2,stage] += atom.IonizRate['rate']
        for stage in range(1,z+1):
            atom = chIons[stage]
#            if type(atom.RecombRate) != types.NoneType:
            evolver[stage,stage] -= atom.RecombRate['rate']
            atom = chIons[stage-1]
            evolver[stage,stage-1] += atom.IonizRate['rate']
#                evolver[stage-1,stage-2] += atom.RecombRate['rate']
        # include normalization of ionization stages
        for stage in range(0,z+1):
            evolver[0,stage] = 1.
        self.Static = evolver
        #
        # -----------------------------------------
        #
    def getEvolver(self,temperature,density,time):
        ''' calculate the evolver matrix for the new temperature'''
        self.NewTemperature = temperature
        z = self.Z
        ionList=[]
        chIons=[]
        evolver = np.zeros((z+1,z+1),'float64')
        for ion in range(1, z+2):
            ions=util.zion2name(z, ion)
            ionList.append(ions)
#            print z, ion, ions
            atom=ch.ion(ions, temperature=temperature, setup=0)
            atom.ionizRate()
            atom.recombRate()
            chIons.append(atom)
        for stage in range(0,z):
#            atom = chIons[stage]
#            if type(atom.IonizRate) != types.NoneType:
            evolver[stage,stage] -= chIons[stage].IonizRate['rate']
#            atom = chIons[stage+1]
            evolver[stage,stage+1] += chIons[stage+1].RecombRate['rate']
#                evolver[stage-2,stage] += atom.IonizRate['rate']
        for stage in range(1,z+1):
#            atom = chIons[stage]
#            if type(atom.RecombRate) != types.NoneType:
            evolver[stage,stage] -= chIons[stage].RecombRate['rate']
#            atom = chIons[stage-1]
            evolver[stage,stage-1] += chIons[stage-1].IonizRate['rate']
#                evolver[stage-1,stage-2] += atom.RecombRate['rate']
        # include normalization of ionization stages
#        for stage in range(0,z+1):
#            evolver[0,stage] = 1.
        # this is the first order Euler solution
#        implicitEvolver = np.invert(np.identity(z+1,'float64') - evolver)
        implicitEvolver = np.dual.inv(np.identity(z+1,'float64') - density*time*evolver)
        self.Evolver = implicitEvolver
#        imvolver = np.asmatrix(np.identity(z+1) - evolver)
#        self.Evolver = np.dual.inv(imvolver)
#        self.Evolver = evolver
        #
        # -----------------------------------------
        #
    def evolve(self,temperature, density, haveEvolver = 0):
        ''' to evolve the plasma to the final temperature'''
        z = self.Z
        self.finalTemperature=temperature
        if not haveEvolver:
            self.getEvolver(temperature,density,1.)
#        diff = np.abs(density*np.dot(self.Evolver,self.Ioneq)/self.Ioneq)
#        diff = np.abs(self.Evolver).max()
#        print ' diffs = ',diff
        maxDiff = 1.
#        print ' max time  = ',1./diff
#        dt = (maxDiff/diff.max())
#        diff = np.abs((np.ones(z+1,'float64') - self.Evolver.diagonal()).max())
        new = np.dot(self.Evolver,self.Ioneq)
        diff = np.where(self.Ioneq > 0.,np.abs(new - self.Ioneq)/self.Ioneq,100.)
        print((' diff = %10.2e'%(diff)))
        dt = maxDiff/diff.min()
        print((' dt = %12.2e'%(dt)))
        self.getEvolver(temperature,density,dt)
        niter = int(10.*maxDiff/diff.min())+1
        print((' niter = %i'%(niter)))
        state = np.zeros((niter+1,z+1),'float64')
        stateTot = np.zeros(niter+1,'float64')
        stateDiff = np.zeros((niter,z+1,z+1),'float64')
        change = np.zeros((niter+1),'float64')
        deltaChange = np.zeros((niter),'float64')
        deltaChange[0] = 1.
        state[0] = self.Ioneq
        stateTot[0] = self.Ioneq.sum()
        previousChange = 1.
        dChange = 10.
        iter = 0
#        for iter in range(niter):
        while (dChange > 1.e-7) and (iter < niter):
            present = state[iter]
#            diff = np.dot(self.Evolver,present)
#            imdiff = np.asmatrix(np.identity(self.Z+1,dtype='float64') - diff)
#            evolver = np.dual.inv(imdiff)
#            diff = diff - diff.sum()
#            new = np.dual.solve(imdiff,present)
            new = np.dot(self.Evolver,present)
            thisdiff = np.where(present > 0.,np.abs(new - present)/present,0.)
#            print iter, change.sum, change
            state[iter+1] = new
            stateTot[iter+1] = new.sum()
#            stateDiff[iter] = imdiff
            change[iter] = thisdiff.sum()
            deltaChange[iter] = np.abs(change[iter] - previousChange)
            dChange = deltaChange[iter]
            previousChange = change[iter]
            iter+=1
        print((' # of iterations = %i'%(iter)))
        self.Temperature = self.finalTemperature
        self.getIoneq()
        state[iter] = self.Ioneq
        self.State = state[:iter].transpose()
        self.StateTot = stateTot[:iter]
        self.StateDiff = stateDiff[:iter]
        self.Change = change[:iter]
        self.DeltaChange = deltaChange[:iter]
        self.Iter = iter







#    def plot(self, ions=None, xr=None, yr=None, oplot=False, label=True, title=True,  bw=False):
#        ''' self.plot(xr=None, yr=None, oplot=False)
#        ions = sequence of ions to be plotted, in spectroscopic notation
#        xr = temperature range, yr = ion fraction range
#        for overplotting:
#        oplot="ioneqfilename" such as 'mazzotta'
#        or oplot=True and a widget will come up so that a file can be selected'''
#        if bw:
#            linestyle=['k-','k--', 'k-.', 'k:']
#        else:
#            linestyle=['b-','r--', 'g-.', 'm:']
#        #
#        if type(ions) == NoneType:
#            ions=range(1, self.Z+2)
#        elif min(ions) < 1 or max(ions) > self.Z+1:
#            ions=range(1, self.Z+2)  #  spectroscopic notation
#        if type(xr) == NoneType:
#            xr=[self.Temperature.min(), self.Temperature.max()]
#        if type(yr) == NoneType:
#            yr=[0.01, 1.1]
#        xyr=list(xr)
#        xyr.extend(list(yr))
#        #
#        iz=ions[0]
#        pl.loglog(self.Temperature, self.Ioneq[iz-1])
#        if label:
#            idx=self.Ioneq[iz-1] == self.Ioneq[iz-1].max()
#            if idx.sum() > 1:
#                jdx=np.arange(len(idx))
#                idx=jdx[idx].max()
#            ann=Ionstage[iz-1]
#            pl.annotate(ann, [self.Temperature[idx], 0.7*self.Ioneq[iz-1, idx]], ha='center')
#        for iz in ions[1:]:
#            pl.plot(self.Temperature, self.Ioneq[iz-1], linestyle[0])
#            if label:
#                idx=self.Ioneq[iz-1] == self.Ioneq[iz-1].max()
#                if idx.sum() > 1:
#                    jdx=np.arange(len(idx))
#                    idx=jdx[idx].mean()
#                ann=Ionstage[iz-1]
#                pl.annotate(ann, [self.Temperature[idx], 0.7*self.Ioneq[iz-1, idx]], ha='center')
#        pl.xlabel('Temperature (K)')
#        pl.ylabel('Ion Fraction')
#        atitle='Chianti Ionization Equilibrium for '+El[self.Z-1].capitalize()
#        #
#        if oplot != False:
#            if type(oplot) == BooleanType:
#                result=self.ioneqRead(ioneqname='',default=False)
#                if result != False:
#                    atitle+='  & '+result['ioneqname'].replace('.ioneq', '')
#                    atitle+=' '+linestyle[0]
#                    for iz in ions:
#                        pl.plot(self.IoneqTemperature, self.IoneqAll[self.Z-1, iz-1],linestyle[0], linestyle[1])
#            elif type(oplot) == StringType:
#                atitle+='  & '+oplot+' '+linestyle[0]
#                atitle+=' '+linestyle[0]
#                result=self.ioneqRead(ioneqname=oplot,default=False)
#                if result != False:
#                    for iz in ions:
#                        pl.plot(self.IoneqTemperature, self.IoneqAll[self.Z-1, iz-1],linestyle[0], linestyle[1])
#            elif type(oplot) == ListType:
#                for iplot in range(len(oplot)):
#                    result=self.ioneqRead(ioneqname=oplot[iplot],default=False)
#                    if result != False:
#                        atitle+='  & '+oplot[iplot]+' '+linestyle[iplot%3]
#                        for iz in ions:
#                            pl.plot(self.IoneqTemperature, self.IoneqAll[self.Z-1, iz-1], linestyle[iplot%4])
#            else:
#                print ' oplot not understood ', oplot
#        if title:
#            pl.title(atitle)
#        pl.axis(xyr)
        #
