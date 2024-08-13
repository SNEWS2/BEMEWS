#!/usr/bin/env python3

import os
import numpy as np

import configparser
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from importlib.resources import files
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

import BEMEWS
import BEMEWS.data

if __name__ == "__main__":
    p = ArgumentParser(description='BEMEWS oscillation probability calculator',
                       formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('-c', '--config', type=str, default='config.ini',
                   help='Config file (INI format)')
    p.add_argument('-s', '--source', type=str,
                   help='Name of supernova source (e.g., Betelguese')
    p.add_argument('-d', '--detector', type=str,
                   help='Detector name (e.g., detector)')
    p.add_argument('-t', '--time', type=str,
                   help='Date and time (ISO format)')
    p.add_argument('--deltam_21', type=float,
                   help='Δm_12 in eV^2')
    p.add_argument('--deltam_32', type=float,
                   help='Δm_32 in eV^2')
    p.add_argument('--theta_12', type=float,
                   help='θ_12 in deg')
    p.add_argument('--theta_13', type=float,
                   help='θ_13 in deg')
    p.add_argument('--theta_23', type=float,
                   help='θ_23 in deg')
    p.add_argument('--delta_cp', type=float,
                   help='δ_CP in deg')
    p.add_argument('--numebins', type=int,
                   help='Number of energy bins')
    p.add_argument('--emin', type=float,
                   help='Minimum energy in MeV')
    p.add_argument('--emax', type=float,
                   help='Maximum energy ax MeV')
    p.add_argument('--accuracy', type=float,
                   help='Accuracy of Runge-Kutta solver')
    p.add_argument('--stepcounterlimit', type=int,
                   help='Output every N steps through the Earth (higher = less output)')

    #- Load arguments from a file, overwriting from the command line if needed
    args = vars(p.parse_args())

    config = configparser.ConfigParser()
    config.read(args['config'])
    defaults = dict()
    for k in config.keys():
        defaults.update(dict(config[k]))
    settings = defaults
    settings.update({k: v for k, v in args.items() if v is not None})

    print(settings)

    # skycoordinates of neutrino source
    source = SkyCoord.from_name(settings['source']) 
    
    # neutrino detector
    detector = EarthLocation.of_site(settings['detector'])

    # time when the supernova occured 
    # the first time option means the neutrinos traveled through the Earth, the second means they did not
    time = Time(settings['time'])
    #time = Time('2021-5-26 14:14:00') - 12*u.hour    

    # altaz of supernovae at detector
    SNaltaz = source.transform_to(AltAz(obstime=time, location=detector)) 

    # class to accumulate input data for calcultion
    ID = BEMEWS.InputDataBEMEWS()

    # assign data fields
    ID.altitude = SNaltaz.alt.deg
    ID.azimuth = SNaltaz.az.deg

    ID.outputfilenamestem = "out/BEMEWS:PREM"

    ID.densityprofile = str(files(BEMEWS.data).joinpath("PREM.rho.dat"))
    ID.electronfraction = str(files(BEMEWS.data).joinpath("PREM.Ye.dat"))

    ID.NE = int(settings['numebins'])
    ID.Emin = float(settings['emin'])
    ID.Emax = float(settings['emax'])

    ID.deltam_21 = float(settings['deltam_21'])
    ID.deltam_32 = float(settings['deltam_32'])
    ID.theta12 = float(settings['theta_12'])
    ID.theta13 = float(settings['theta_13'])
    ID.theta23 = float(settings['theta_23'])
    ID.deltaCP = float(settings['delta_cp'])

    ID.accuracy = float(settings['accuracy'])

    # if set to True the BEMEWS module will output files in the 'out' directory
    # The stepcounterlimit controls how often output is written. The larger the number, the less often it happens. 
    ID.outputflag = settings['output'] in ['True', 'true']
    ID.stepcounterlimit = int(settings['stepcounterlimit'])
    
    if ID.outputflag:
        os.makedirs('out', exist_ok=True)

    # do the calculation. The return is a four dimensional array of transition probabilities nu_alpha -> nu_i: 
    # index order is matter/antimatter, energy, i, alpha
    Pfm = BEMEWS.Run(ID)

    print("finished")
