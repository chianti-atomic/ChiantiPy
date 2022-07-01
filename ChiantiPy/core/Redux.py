import os.path
from ChiantiPy.base import ionTrails
from ChiantiPy.base import specTrails


class redux(ionTrails, specTrails):
    """ a class for restoring save multi-ion saveData files
    """
    def __init__(self,  filename,  verbose=False):
        """

        :param filename: filename of a file saved by the saveData method
        :type filename: str

        :param verbose: if true, prints additional informations
        :type verbose:  bool

        the redux class inherets the following useful methods:
            intensityList
            intensityPlot
            intensityRatio
            intensityRatioSave

            convolve
            lineSpectrumPlot
            restoreData
            spectrumPlot

        """
        if not os.path.isfile(filename):
            print(' print restore must be set to a valid filename')
            return
        elif os.path.isfile(filename):
            self.restoreData(filename)

        if verbose:
            for aname in self.__dict__.keys():
                print(' restored attribute %s'%(aname))
