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

BEMEWS is a python module for calculating the Earth-Matter Effect (EME) upon neutrinos. It is a standalone module but it can also be used with SNEWPY 2.0 to add EMEs to supernova neutrino signals. 

# Statement of need

If the neutrinos from a source pass through the Earth before reaching a detector then an imprint from the matter can be left on the signal. This imprint depends upon the profle of the matter along the trajectory between source and detector. BEMEWS (Better Earth Matter Effects With SNEWPY) is a python module that calculates the Earth-Matter Effect for a given source and detector location, neutrino energy and mixing parameters. The EME is determined by solving the Schrodinger equation for the neutrino evolution operator using a Hamiltonian composed of the 'vacuum term' and the MSW 'matter term'. Further details of the calculation can be found in the appendix of the in the SNEWPY 2.0 documentation. The output is in the form of tables of the probability that a neutrino of a given initial state nu_j is detected as state nu_i. Output files of the probabilities as a function of distance traveled along the trajectroy are generated for each neutrino energy, and of the probabilities as a function of neutrino energy at the detector location. 

# How to use BEMEWS

The BEMEWS.py example script in the doc folder shows how to use the module in standalone mode. BEMEWS can also be imported into the SNEWPY 2.0  software and used as a `modifier' in a TransformationChain neutrino flavor transformation prescription. How to include BEMEWS with SNEWPY will be documented in the SNEWPY documentation. Note that BEMEWS uses the PREM as the density profile for the Earth: this can be changed by altering the contents of the file in the src/data folder. 

As shown in the BEMEWS.py example script, to use BEMEWS the user must provide the altitude and azimuth of the source relative to the detector. For supernova neutrinos from a particular sky location, these can be determined from an ASTROPY AltAz object. To facilitate computing the AltAz objects for specific sky and detector locations, SNEWS has worked with ASTROPY so that ASTROPY now includes many neutrino detectors in its list of EarthLocation classes that can be referenced by name. At the present time these detectors are: IceCube, NOvA, HALO, SNO+, ORCA, ARCA, HyperK and SuperK. In the BEMEWS.py example script, an ASTROPY AltAz object is made for the combination of the explosion of Betelgeuse on May 26th 2021 at 14:14:00 UT, and the SuperK detector.  

# Acknowledgements

This work is supported by the National Science Foundation “Windows on the Universe: the Era of Multi-Messenger Astrophysics” Program: “WoU-MMA:
Collaborative Research: A Next-Generation SuperNova Early Warning System for Multimessenger Astronomy” through Grant Nos. 1914448, 1914409, 1914447,
1914418, 1914410, 1914416, and 1914426. This work is also supported at NC State by U.S. Department of Energy grant DE-FG02-02ER41216, and at King’s College London by STFC.

# References

