'''
Module for collecting various top-level Chianti data

Descriptions for `keywordArgs`:

- `temperature` : temperature (in K)
- `eDensity` : electron density (in :math:`\mathrm{cm}^{-3}`)
- `hDensity` : hydrogen density (in :math:`\mathrm{cm}^{-3}`)
- `pDensity` : proton density (in :math:`\mathrm{cm}^{-3}`)
- `radTemperature` : radiation temperature of central source (in K)
- `rStar` : distance of the plasma from the source (in units of the source's radius)
- `distance` : distance from the central source
'''
import os

import ChiantiPy.tools.io as chio

###
xuvtop = os.environ['XUVTOP']
#chInteractive=1
Defaults = chio.defaultsRead()
Ip = chio.ipRead()
MasterList = chio.masterListRead()
IoneqAll = chio.ioneqRead(ioneqname = Defaults['ioneqfile'])
# gets the version of the CHIANTI database
ChiantiVersion = chio.versionRead()
keywordArgs = ['temperature','eDensity','hDensity', 'pDensity','radTemperature', 'rStar', 'distance']
#
AbundanceDefault = chio.abundanceRead(abundancename = Defaults['abundfile'])
abunddir = os.path.join(xuvtop,'abundance')
filelist = os.listdir(abunddir)
#
AbundanceList = []
for one in filelist:
    fname = os.path.join(abunddir,one)
    if os.path.isfile(fname):
        AbundanceList.append(os.path.splitext(one)[0])
#for one in abundList:
#    print(one)
Abundance = {AbundanceList[0]:chio.abundanceRead(abundancename = AbundanceList[0])}
for one in AbundanceList[1:]:
    Abundance[one] = chio.abundanceRead(abundancename = one)
