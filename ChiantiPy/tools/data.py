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

Parameters
----------

XUVTOP : str
    the root of the CHIANTI database
Defaults : dict
    default values used by ChiantiPy, can also be set by an optional HOME/.chianti/chiantirc file
Ip : np.ndarray
    the ionization potentials for all ionization stages up to Zn, in eV
MasterList : list
    the CHIANTI style names of all ions in the CHIANTI database
IoneqAll : dict
    a dict containing the ionization equilibrium values for the default ionization file
ChiantiVersion : str
    the version of the CHIANTI database
AbundanceDefault : dict
    the elemental abundances in the default abundance file
AbundanceList : list
    the names of all abundance files included in the CHIANTI database
GrndLevels : list
    the number of levels that should be considered in an ionization calculation
'''
import os
import glob
import traceback
import warnings

import ChiantiPy.tools.io as chio

try:
    Xuvtop = os.environ['XUVTOP']
    Defaults = chio.defaultsRead()
    Ip = chio.ipRead()
    MasterList = chio.masterListRead()
    IoneqAll = chio.ioneqRead(ioneqName=Defaults['ioneqfile'])
    ChiantiVersion = chio.versionRead()
    keywordArgs = ['temperature', 'eDensity', 'hDensity', 'pDensity', 'radTemperature',
                   'rStar', 'distance']
    AbundanceDefault = chio.abundanceRead(abundancename=Defaults['abundfile'])

    AbundanceList = []
    for fname in glob.glob(os.path.join(Xuvtop, 'abundance', '*.abund')):
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
            'Cannot find the CHIANTI atomic database at {}. You will not be able to access any data from the CHIANTI database.'.format(Xuvtop))
