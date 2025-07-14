import unittest

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

from importlib import import_module
from importlib.resources import files
import numpy as np
import sys
import os

import BEMEWS

class TestInit(unittest.TestCase):

    def test_run_no_output(self):
        """Test run with no output data.
        """
        _id = BEMEWS.InputDataBEMEWS()

        source = SkyCoord.from_name('Betelgeuse')
        detector = EarthLocation.of_site('SuperK')
        time = Time('2021-05-26 14:14:00')
        snaltaz = source.transform_to(AltAz(obstime=time, location=detector))

        _id.altitude = snaltaz.alt.deg
        _id.azimuth = snaltaz.az.deg
        _id.outputfilenamestem = 'out/BEMEWS:PREM'

        _id.densityprofile = str(files(BEMEWS.data).joinpath('PREM.rho.dat'))
        _id.electronfraction = str(files(BEMEWS.data).joinpath('PREM.Ye.dat'))

        _id.deltam_21 = 7.69e-5
        _id.deltam_32 = 2.43e-3
        _id.theta12 = 34.4
        _id.theta13 = 9
        _id.theta23 = 45
        _id.deltaCP = 0

        _id.NE = 2
        _id.Emin = 1
        _id.Emax = 10

        _id.outputflag = False
        _id.ecsvformat = False
        _id.accuracy = 1.01e-9
        _id.stepcounterlimit = 1

        Pfm = np.asarray(BEMEWS.Run(_id))

        self.assertTrue(
            np.all(
                np.isclose(Pfm[0,0,:],
                    np.asarray([[0.66414735, 0.31137679, 0.02447587],
                                [0.24085689, 0.27138279, 0.48776032],
                                [0.09499576, 0.41724042, 0.48776382]]))))

        self.assertTrue(
            np.all(
                np.isclose(Pfm[1,0,:],
                    np.asarray([[0.66533681, 0.31019902, 0.02446416],
                                [0.24016784, 0.27205631, 0.48777586],
                                [0.09449535, 0.41774467, 0.48775998]]))))
