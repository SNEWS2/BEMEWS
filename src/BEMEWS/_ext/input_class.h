
#include "BEMEWS.h"

#ifndef input_class_H
#define input_class_H

// **************************************

struct InputDataBEMEWS;

void Profile_loader(InputDataBEMEWS ID,std::string &outputfilenamestem);
void Neutrino_loader(InputDataBEMEWS ID,std::string &outputfilenamestem);

// **********************************************************
// **********************************************************
// **********************************************************

struct InputDataBEMEWS 
       { double altitude, azimuth; // in decimal degrees, altitude is negative for angles below the horizon
         std::string outputfilenamestem;
         std::string densityprofile;
         std::string electronfraction;

         int NE; // number of energies
         double Emin, Emax; // in MeV
         double deltam_21, deltam_32; // in eV^2	
         double theta12, theta13, theta23, deltaCP; // all in degrees
         double accuracy;
         double stepcounterlimit; // how often it outputs data
         bool outputflag; // whether the code outputs data as it does the integration
         bool ecsvformat; // whether the output is in ECSV or tab-separated text files

         InputDataBEMEWS(void) {;}  
        };

#endif
