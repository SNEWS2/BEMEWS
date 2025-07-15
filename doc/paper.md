---
title: 'BEMEWS: Better Earth Matter Effects With SNEWPY'
tags:
  - Python
  - astronomy
  - supernova
  - neutrinos
authors:
  - name: Segev BenZvi
    orcid: 0000-0001-5537-4710
    affiliation: 1
  - name: Marta Colomer Molla
    orcid: 0000-0003-1801-8121
    affiliation: 2
  - name: Alec Habig
    affiliation: 3
    orcid: 0000-0002-1018-9383
  - name: James P. Kneller^[Corresponding author]
    orcid: 0000-0002-3502-3830
    affiliation: 4
  - name: Jost Migenda
    orcid: 0000-0002-5350-8049
    affiliation: 5
  - name: Kate Scholberg
    orcid: 0000-0002-7007-2021
    affiliation: 6
  - name: Andrey Sheshukov
    affiliation: 7
    orcid: 0000-0001-5128-9279
  - name: Jeff Tseng
    affiliation: 8
    orcid: 0000-0003-1731-5853
affiliations:
  - name: University of Rochester, Rochester, NY, USA
    index: 1
  - name: Université Libre de Bruxelles, Brussels, Belgium
    index: 2
  - name: University of Minnesota Duluth, Duluth, MN, USA
    index: 3
  - name: NC State University, Raleigh, NC, USA
    index: 4    
  - name: King’s College London, London, UK
    index: 5
  - name: Duke University, Durham, NC, USA
    index: 6
  - name: Joint Institute for Nuclear Research, Dubna, Russia
    index: 7
  - name: Oxford University, Oxford, UK
    index: 8
date: 1 September 2025
bibliography: paper.bib

---

# Summary

BEMEWS is a python module for calculating the Earth-matter effect on neutrino flavor transformations. It is a standalone module that SNEWPY uses to compute the Earth-Matter Effect for supernova neutrinos. The BEMEWS_example.py script shows how to use the module in standalone mode. The EarthMatter flavor transformation class in SNEWPY is essentially the same script but with the output options turned off.


# Statement of need

If the neutrinos from a supernova pass through the Earth before reaching a detector then an imprint can be left on the signal. This imprint depends upon the location of the supernova on the sky relative to the detector. BEMEWS (Better Earth Matter Effects With SNEWPY) is a python module that calculates the Earth Matter Effect for a given sky location of the supernova, Earth location for the detector, neutrino energy and mixing paramaters. It can be run as a standalone code or it can imported into the SNEWPY software and used as a `modifier' in a TransformationChain neutrino flavor transformation prescription. 

The package, written in Python, is built upon NUMPY [@harris2020array] and SCIPY [@Virtanen:2019joe], and makes use of  ASTROPY [@Astropy:2013muo;
@Price-Whelan:2018hus] for angle conversions, sky location of well-known progenitors, and Earth locations of neutrinos detecors.

SNEWPY will function without access to the BEMEWS module but of course the EarthMatter modifier will not be available. Note that BEMEWS 
uses the PREM as the density profile for the Earth: this can be changed by altering the contents of the datafile. 

To use the module the user must create an ASTROPY AltAz object with the altitude-azimuth (AltAz) of the supernova at the detector. To facilitate computing the AltAz objects, ASTROPY now includes many neutrino detectors in its list of EarthLocation classes that can be referenced by name. At the present time these detectors are: IceCube, NOvA, HALO, SNO+, ORCA, ARCA, HyperK and SuperK. 
Once made, the AstroPY AltAz object can be input into the EarthMatter prescription which itself is then input into the TransformationChain prescription. An example script is provided in which an ASTROPY AltAz object for the SuperK detector if Betelgeuse exploded on May 26th 2021 at 23:14:00 local time is created. 

# Acknowledgements

This work is supported by the National Science Foundation “Windows on the Universe: the Era of Multi-Messenger Astrophysics” Program: “WoU-MMA:
Collaborative Research: A Next-Generation SuperNova Early Warning System for Multimessenger Astronomy” through Grant Nos. 1914448, 1914409, 1914447,
1914418, 1914410, 1914416, and 1914426. This work is also supported at NC State by U.S. Department of Energy grant DE-FG02-02ER41216, and at King’s College London by STFC.

# References

