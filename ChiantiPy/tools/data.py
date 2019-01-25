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
import glob
import traceback
import warnings

import ChiantiPy.tools.io as chio

try:
    xuvtop = os.environ['XUVTOP']
    Defaults = chio.defaultsRead()
    Ip = chio.ipRead()
    MasterList = chio.masterListRead()
    IoneqAll = chio.ioneqRead(ioneqName=Defaults['ioneqfile'])
    ChiantiVersion = chio.versionRead()
    keywordArgs = ['temperature', 'eDensity', 'hDensity', 'pDensity', 'radTemperature',
                   'rStar', 'distance']
    AbundanceDefault = chio.abundanceRead(abundancename=Defaults['abundfile'])

    AbundanceList = []
    for fname in glob.glob(os.path.join(xuvtop, 'abundance', '*.abund')):
        AbundanceList.append(os.path.splitext(os.path.basename(fname))[0])
    Abundance = {abundance: chio.abundanceRead(abundancename=abundance)
                 for abundance in AbundanceList}
    GrndLevels = chio.grndLevelsRead()
except (KeyError, IOError) as e:
    print(traceback.format_exc())
    if isinstance(e, KeyError):
        warnings.warn(
            'XUVTOP environment variable not set. You will not be able to access any data from the CHIANTI database.')
    else:
        warnings.warn(
            'Cannot find the CHIANTI atomic database at {}. You will not be able to access any data from the CHIANTI database.'.format(xuvtop))
