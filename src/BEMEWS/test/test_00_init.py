import unittest

from importlib import import_module
import numpy as np
import sys

import BEMEWS
import BEMEWS.data

class TestInit(unittest.TestCase):
    def test_import(self):
        """Test import.
        """
        import_module('BEMEWS')
        self.assertTrue('BEMEWS' in sys.modules)

    def test_inputdata(self):
        """Test initialization and access of the BEMEWS settings struct.
        """
        _id = BEMEWS.InputDataBEMEWS()

        _id.altitude = 10.1
        self.assertTrue(np.isclose(_id.altitude, 10.1))

        _id.azimuth = 99.9
        self.assertTrue(np.isclose(_id.azimuth, 99.9))

        _id.outputfilenamestem = 'out/BEMEWS:PREM'
        self.assertTrue(_id.outputfilenamestem == 'out/BEMEWS:PREM')

        _id.deltam_21 = 7.69e-5
        self.assertTrue(np.isclose(_id.deltam_21, 7.69e-5))

        _id.deltam_32 = 2.43e-3
        self.assertTrue(np.isclose(_id.deltam_32, 2.43e-3))

        _id.theta12 = 34.4
        self.assertTrue(np.isclose(_id.theta12, 34.4))

        _id.theta13 = 9
        self.assertTrue(np.isclose(_id.theta13, 9))

        _id.theta23 = 45
        self.assertTrue(_id.theta23 == 45)

        _id.deltaCP = 0
        self.assertTrue(_id.deltaCP == 0)

        _id.NE = 296
        self.assertTrue(_id.NE == 296)

        _id.Emin = 1
        self.assertTrue(_id.Emin == 1)

        _id.Emax = 60
        self.assertTrue(_id.Emax == 60)

        _id.outputflag = False
        self.assertTrue(not _id.outputflag)

        _id.ecsvformat = False
        self.assertTrue(not _id.ecsvformat)

        _id.accuracy = 1.01e-9
        self.assertTrue(np.isclose(_id.accuracy, 1.01e-9))

        _id.stepcounterlimit = 1
        self.assertTrue(_id.stepcounterlimit == 1)
