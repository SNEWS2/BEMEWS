
#include "update.h"

// *****************************************************************************

using std::complex;
using std::array;
using std::vector;

// ***************************** UpdateSm ***************************************

vector<vector<MATRIX<complex<double>,NF,NF> > > UpdateSm(double lambdaminus,double lambdaplus,vector<vector<array<double,NY> > > &Y,vector<vector<MATRIX<complex<double>,NF,NF> > > &Smprior)
          { array<MATRIX<complex<double>,NF,NF>,NM> VfMSW;
            MATRIX<complex<double>,NF,NF> Hf,Hfbar;
            MATRIX<complex<double>,NF,NF> UU,UUbar; 

            array<double,NF> kk,kkbar,dkk,dkkbar;

            vector<vector<MATRIX<complex<double>,NF,NF> > > Sm(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));

            double r, rrho, YYe;

            // multiply SMSW by the mixing matrix at rminus
            r=sqrt(RE*RE+lambdaminus*lambdaminus+2.*RE*lambdaminus*sin(altitude));
            rrho = rho(r);
            YYe = Ye(r); 
    
            VfMSW[nu][e][e] = Ve(rrho,YYe); 
            VfMSW[nu][mu][mu] = Vmu(rrho,YYe); 
            VfMSW[nu][tau][tau] = Vtau(rrho,YYe); 

            VfMSW[antinu] = -VfMSW[nu];

            int i;
	    #pragma omp parallel for schedule(static) private(Hf,Hfbar,kk,kkbar,dkk,dkkbar,UU,UUbar)
	    for(i=0;i<=NE-1;i++)
               { Hf=HfV[nu][i]+VfMSW[nu];
                 kk=k(Hf);
                 dkk=deltak(kk);
                 UU=MixingMatrix(Hf,kk,dkk);
            	 Sm[nu][i]=UU * W(Y[nu][i]) * B(Y[nu][i]) * Smprior[nu][i];

            	 Hfbar=HfV[antinu][i]+VfMSW[antinu];
            	 kkbar=kbar(Hfbar);
            	 dkkbar=deltakbar(kkbar);
                 UUbar=MixingMatrix(Hfbar,kkbar,dkkbar);
            	 Sm[antinu][i]=UUbar * W(Y[antinu][i]) * B(Y[antinu][i]) * Smprior[antinu][i];
		}

            // multiply SMSW by the adjoint of the mixing matrix at rplus
            r=sqrt(RE*RE+lambdaplus*lambdaplus+2.*RE*lambdaplus*sin(altitude));
            rrho = rho(r);
            YYe = Ye(r); 
    
            VfMSW[nu][e][e] = Ve(rrho,YYe); 
            VfMSW[nu][mu][mu] = Vmu(rrho,YYe); 
            VfMSW[nu][tau][tau] = Vtau(rrho,YYe); 

            VfMSW[antinu] = -VfMSW[nu];

	    #pragma omp parallel for schedule(static) private(Hf,Hfbar,kk,kkbar,dkk,dkkbar,UU,UUbar)
	    for(i=0;i<=NE-1;i++)
               { Hf=HfV[nu][i]+VfMSW[nu];
                 kk=k(Hf);
             	 dkk=deltak(kk);
                 UU=MixingMatrix(Hf,kk,dkk);  
            	 Sm[nu][i]=Adjoint(UU)*MATRIX<complex<double>,NF,NF>(Sm[nu][i]); 

            	 Hfbar=HfV[antinu][i]+VfMSW[antinu];
                 kkbar=kbar(Hfbar);
            	 dkkbar=deltakbar(kkbar);
                 UUbar=MixingMatrix(Hfbar,kkbar,dkkbar);
            	 Sm[antinu][i]=Adjoint(UUbar)*MATRIX<complex<double>,NF,NF>(Sm[antinu][i]);
	        }

            return Sm;
        }





