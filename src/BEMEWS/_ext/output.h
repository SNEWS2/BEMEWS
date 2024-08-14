
#include "BEMEWS.h"

struct InputDataBEMEWS;

#ifndef output_H
#define output_H

extern std::vector<std::string> fPvsrfilename;

void Initialize_Output(std::string outputfilenamestem, std::ofstream &fPvslambda, std::ofstream &fHvslambda, const InputDataBEMEWS &id);

void Close_Output(std::ofstream &fHvslambda);

void Output_Pvslambda(bool firsttime, bool lasttime, std::ofstream &fPvslambda, double lambda, std::vector<std::vector<std::array<double, NY> > > &Y, std::vector<std::vector<MATRIX<std::complex<double>, NF, NF> > > &Scumulative, const InputDataBEMEWS &id);

void Output_PvsE(bool lasttime, std::ofstream &fPvsE, std::string outputfilenamestem, double lambda, std::vector<std::vector<std::array<double, NY> > > &Y, std::vector<std::vector<MATRIX<std::complex<double>, NF, NF> > > &Scumulative, const InputDataBEMEWS &id);

void Output_Hvslambda(bool firsttime, bool lasttime, std::ofstream &fHvslambda, double lambda, std::vector<std::vector<std::array<double, NY> > > &Y, std::vector<std::vector<MATRIX<std::complex<double>, NF, NF> > > &Scumulative, const InputDataBEMEWS &id);

#endif
