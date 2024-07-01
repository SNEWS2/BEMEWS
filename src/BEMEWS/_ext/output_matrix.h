
#include "BEMEWS.h"

#ifndef output_matrix_H
#define output_matrix_H

void Pfm(double lambda,std::vector<std::vector<std::array<double,NY> > > &Y,std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > &Scumulative,std::vector<std::vector<std::vector<std::vector<double> > > > &PPfm);

void Pmf(double r,std::vector<std::vector<std::array<double,NY> > > &Y,std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > &Scumulative,std::vector<std::vector<std::vector<std::vector<double> > > > &PPmf);

#endif
