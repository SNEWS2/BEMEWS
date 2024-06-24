
#include "EMEWS.h"

#ifndef update_H
#define update_H

std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > UpdateSm(double lambdaminus,double lambdaplus,std::vector<std::vector<std::array<double,NY> > > &Y,std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > &Smprior);

#endif
