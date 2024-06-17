
#include "EMEWS.h"

#ifndef output_H
#define output_H

extern std::vector<std::string> fPvsrfilename;

void Initialize_Output(std::string outputfilenamestem,std::ofstream &fPvslambda,std::ofstream &fHvslambda);

void Close_Output(std::ofstream &fHvslambda);

void Output_Pvslambda(bool firsttime,std::ofstream &fPvslambda,double lambda,std::vector<std::vector<std::array<double,NY> > > &Y,std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > &Scumulative);

void Output_PvsE(std::ofstream &fPvsE,std::string outputfilenamestem,double lambda,std::vector<std::vector<std::array<double,NY> > > &Y,std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > &Scumulative);

void Output_Hvslambda(bool firsttime,std::ofstream &fHvslambda,double lambda,std::vector<std::vector<std::array<double,NY> > > &Y,std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > &Scumulative);

#endif
