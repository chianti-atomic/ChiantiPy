import numpy as np
import ChiantiPy.tools.io as io
from scipy.interpolate import splev, splrep

class ioneqOne(object):
    """
    Base class for `ChiantiPy.core.ion` and `ChiantiPy.core.continuum`
    """
    def ioneqOne(self):
        '''
        Provide the ionization equilibrium for the selected ion as a function of temperature.

        returned in self.IoneqOne
        '''
        #
        if hasattr(self, 'Temperature'):
            temperature = self.Temperature
        else:
            return
        #
        if hasattr(self, 'IoneqAll'):
            ioneqAll = self.IoneqAll
        else:
            ioneqAll = io.ioneqRead(ioneqname = self.Defaults['ioneqfile'])
            self.ioneqAll = self.IoneqAll
        #
        ioneqTemperature = ioneqAll['ioneqTemperature']
        Z = self.Z
        stage = self.Ion
        Dielectronic = self.Dielectronic
        ioneqOne = np.zeros_like(temperature)
        #
        thisIoneq = ioneqAll['ioneqAll'][Z-1,stage-1 + Dielectronic].squeeze()
        gioneq = thisIoneq > 0.
        goodt1 = self.Temperature >= ioneqTemperature[gioneq].min()
        goodt2 = self.Temperature <= ioneqTemperature[gioneq].max()
        goodt = np.logical_and(goodt1,goodt2)
        y2 = splrep(np.log(ioneqTemperature[gioneq]),np.log(thisIoneq[gioneq]),s=0)
        #
        if goodt.sum() > 0:
            if self.Temperature.size > 1:
                gIoneq = splev(np.log(self.Temperature[goodt]),y2)   #,der=0)
                ioneqOne[goodt] = np.exp(gIoneq)
            else:
                gIoneq = splev(np.log(self.Temperature),y2)
                ioneqOne = np.exp(gIoneq)*np.ones(self.NTempDens, np.float64)
            self.IoneqOne = ioneqOne
        else:
            self.IoneqOne = np.zeros_like(self.Temperature)
