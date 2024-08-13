
#include "output.h"

// *******************************************************************************

using std::string;
using std::stringstream;

using std::ofstream;

using std::complex;
using std::array;
using std::vector;

using namespace prefixes;

// ********************************************************************************

std::vector<std::string> fPvslambdafilename;

void Initialize_Output(string outputfilenamestem,ofstream &fPvslambda,ofstream &fHvslambda, bool ecsvformat)
         { stringstream filename; 

           fPvslambdafilename=vector<string>(NE);

           for(int i=0;i<=NE-1;i++)
              { filename.str("");
                if(NE>1){ filename << outputfilenamestem << string(":E=") << ((NE-1.-i)*EminMeV+i*EmaxMeV)/(NE-1.) << string("MeV:Pvslambda.dat");}
                else{ filename << outputfilenamestem << string(":E=") << EminMeV << string("MeV:Pvslambda.dat");}
                fPvslambdafilename[i]=filename.str();
                fPvslambda.open(fPvslambdafilename[i].c_str()); fPvslambda.close(); // clears the file
               }

           filename.str("");
           filename << outputfilenamestem <<string(":Hvslambda.dat");
           fHvslambda.open((filename.str()).c_str());
           fHvslambda.precision(12);
          }

// ******************************************************

void Close_Output(ofstream &fHvslambda)
         { fHvslambda.close();}

// ******************************************************

void Output_Pvslambda(bool firsttime,bool lasttime,ofstream &fPvslambda,double lambda,vector<vector<array<double,NY> > > &Y,vector<vector<MATRIX<complex<double>,NF,NF> > > &Scumulative, bool ecsvformat)
      { array<MATRIX<complex<double>,NF,NF>,NM> VfMSW, dVfMSWdlambda;

        double r = sqrt( RE*RE + lambda*lambda - 2.*RE*lambda*sin(-altitude) );
        if(r > RE){ r = RE;}       
        double rrho=rho(r);
        double YYe=Ye(r);

        VfMSW[nu][e][e]=Ve(rrho,YYe);
        VfMSW[nu][mu][mu]=Vmu(rrho,YYe);
        VfMSW[nu][tau][tau]=Vtau(rrho,YYe);
        VfMSW[antinu]=-VfMSW[nu];

        vector<vector<MATRIX<complex<double>,NF,NF> > > Hf(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));
        vector<vector<MATRIX<complex<double>,NF,NF> > > UU(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));        
        vector<vector<MATRIX<complex<double>,NF,NF> > > Sa(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)), Sm(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)), Smf(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)), Sfm(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)), Sf(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)); 

        vector<vector<array<double,NF> > > kk(NM,vector<array<double,NF> >(NE));
        vector<vector<array<double,NF> > > dkk(NM,vector<array<double,NF> >(NE));
      
        int i;
        #pragma omp parallel for schedule(static)
        for(i=0;i<=NE-1;i++)
           { Hf[nu][i]=HfV[nu][i] + VfMSW[nu];
             kk[nu][i]=k(Hf[nu][i]);
                dkk[nu][i]=deltak(kk[nu][i]);
             UU[nu][i]=MixingMatrix(Hf[nu][i],kk[nu][i],dkk[nu][i]);

             Sa[nu][i] = W(Y[nu][i]) * B(Y[nu][i]);                

             if(firsttime==false && lasttime==false){ Sm[nu][i] = Sa[nu][i] * Scumulative[nu][i];}
             else{ Sm[nu][i] = Adjoint(UV[nu])*UU[nu][i] * Sa[nu][i] * Scumulative[nu][i];}
                          
             Smf[nu][i]= Sm[nu][i] * Adjoint(UV[nu]);
             
             if(lasttime==false){ 
                Sfm[nu][i] = UU[nu][i] * Sm[nu][i];
                Sf[nu][i] = UU[nu][i] * Smf[nu][i];
               } 
             else{ Sfm[nu][i] = UV[nu] * Sm[nu][i];
                   Sf[nu][i] = UV[nu] * Smf[nu][i];
                  }
             // *******

             Hf[antinu][i]=HfV[antinu][i] + VfMSW[antinu];
                 kk[antinu][i]=kbar(Hf[antinu][i]);
             dkk[antinu][i]=deltakbar(kk[antinu][i]);
             UU[antinu][i]=MixingMatrix(Hf[antinu][i],kk[antinu][i],dkk[antinu][i]);
      
             Sa[antinu][i] = W(Y[antinu][i]) * B(Y[antinu][i]);

             if(firsttime==false && lasttime==false){ Sm[antinu][i] = Sa[antinu][i] * Scumulative[antinu][i];}
             else{ Sm[antinu][i] = Adjoint(UV[antinu])*UU[antinu][i] * Sa[antinu][i] * Scumulative[antinu][i];}
             
             Smf[antinu][i]= Sm[antinu][i] * Adjoint(UV[antinu]);
             
             if(lasttime==false)
               { Sfm[antinu][i] = UU[antinu][i] * Sm[antinu][i];
                     Sf[antinu][i] = UU[antinu][i] * Smf[antinu][i];
                    }
             else{ Sfm[antinu][i] = UV[antinu] * Sm[antinu][i];
                       Sf[antinu][i] = UV[antinu] * Smf[antinu][i];
                      }         
            }

        for(i=0;i<=NE-1;i++)
           { fPvslambda.open(fPvslambdafilename[i].c_str(),std::ofstream::app);
             fPvslambda.precision(12);
             
             if(firsttime==true){
                fPvslambda<<"lambda [cm] \t r [cm]";

                fPvslambda<<"\t P11 \t P12 \t P13 \t P21 \t P22 \t P23 \t P31 \t P32 \t P33";
                fPvslambda<<"\t Pbar11 \t Pbar12 \t Pbar13 \t Pbar21 \t Pbar22 \t Pbar23 \t Pbar31 \t Pbar32 \t Pbar33";

                fPvslambda<<"\t Pe1 \t Pe2 \t Pe3 \t Pmu1 \t Pmu2 \t Pmu3 \t Ptau1 \t Ptau2 \t Ptau3";
                fPvslambda<<"\t Pbare1 \t Pbare2 \t Pbare3 \t Pbarmu1 \t Pbarmu2 \t Pbarmu3 \t Pbartau1 \t Pbartau2 \t Pbartau3";

                fPvslambda<<"\t Pee \t Pemu \t Petau \t Pmue \t Pmumu \t Pmutau \t Ptaue \t Ptaumu \t Ptautau";
                fPvslambda<<"\t Pbaree \t Pbaremu \t Pbaretau \t Pbarmue \t Pbarmumu \t Pbarmutau \t Pbartaue \t Pbartaumu \t Pbartautau";
               }             

             if(firsttime==true){ fPvslambda<<"\n"<<lambdamin0<<"\t"<<r;}
             if(lasttime==true){ fPvslambda<<"\n"<<lambdamax0<<"\t"<<r;}
             if(firsttime==false && lasttime==false){ fPvslambda<<"\n"<<lambda<<"\t"<<r;}

             fPvslambda<<"\t"<<norm(Sm[nu][i][0][0])<<"\t"<<norm(Sm[nu][i][0][1])<<"\t"<<norm(Sm[nu][i][0][2]);
             fPvslambda<<"\t"<<norm(Sm[nu][i][1][0])<<"\t"<<norm(Sm[nu][i][1][1])<<"\t"<<norm(Sm[nu][i][1][2]);
             fPvslambda<<"\t"<<norm(Sm[nu][i][2][0])<<"\t"<<norm(Sm[nu][i][2][1])<<"\t"<<norm(Sm[nu][i][2][2]);

             fPvslambda<<"\t"<<norm(Sm[antinu][i][0][0])<<"\t"<<norm(Sm[antinu][i][0][1])<<"\t"<<norm(Sm[antinu][i][0][2]);
             fPvslambda<<"\t"<<norm(Sm[antinu][i][1][0])<<"\t"<<norm(Sm[antinu][i][1][1])<<"\t"<<norm(Sm[antinu][i][1][2]);
             fPvslambda<<"\t"<<norm(Sm[antinu][i][2][0])<<"\t"<<norm(Sm[antinu][i][2][1])<<"\t"<<norm(Sm[antinu][i][2][2]);

             fPvslambda<<"\t"<<norm(Sfm[nu][i][e][0])<<"\t"<<norm(Sfm[nu][i][e][1])<<"\t"<<norm(Sfm[nu][i][e][2]);
             fPvslambda<<"\t"<<norm(Sfm[nu][i][mu][0])<<"\t"<<norm(Sfm[nu][i][mu][1])<<"\t"<<norm(Sfm[nu][i][mu][2]);
             fPvslambda<<"\t"<<norm(Sfm[nu][i][tau][0])<<"\t"<<norm(Sfm[nu][i][tau][1])<<"\t"<<norm(Sfm[nu][i][tau][2]);

             fPvslambda<<"\t"<<norm(Sfm[antinu][i][e][0])<<"\t"<<norm(Sfm[antinu][i][e][1])<<"\t"<<norm(Sfm[antinu][i][e][2]);
             fPvslambda<<"\t"<<norm(Sfm[antinu][i][mu][0])<<"\t"<<norm(Sfm[antinu][i][mu][1])<<"\t"<<norm(Sfm[antinu][i][mu][2]);
             fPvslambda<<"\t"<<norm(Sfm[antinu][i][tau][0])<<"\t"<<norm(Sfm[antinu][i][tau][1])<<"\t"<<norm(Sfm[antinu][i][tau][2]);

             fPvslambda<<"\t"<<norm(Sf[nu][i][e][e])<<"\t"<<norm(Sf[nu][i][e][mu])<<"\t"<<norm(Sf[nu][i][e][tau]);
             fPvslambda<<"\t"<<norm(Sf[nu][i][mu][e])<<"\t"<<norm(Sf[nu][i][mu][mu])<<"\t"<<norm(Sf[nu][i][mu][tau]);
             fPvslambda<<"\t"<<norm(Sf[nu][i][tau][e])<<"\t"<<norm(Sf[nu][i][tau][mu])<<"\t"<<norm(Sf[nu][i][tau][tau]);

             fPvslambda<<"\t"<<norm(Sf[antinu][i][e][e])<<"\t"<<norm(Sf[antinu][i][e][mu])<<"\t"<<norm(Sf[antinu][i][e][tau]);
             fPvslambda<<"\t"<<norm(Sf[antinu][i][mu][e])<<"\t"<<norm(Sf[antinu][i][mu][mu])<<"\t"<<norm(Sf[antinu][i][mu][tau]);
             fPvslambda<<"\t"<<norm(Sf[antinu][i][tau][e])<<"\t"<<norm(Sf[antinu][i][tau][mu])<<"\t"<<norm(Sf[antinu][i][tau][tau]);

             fPvslambda.flush();
             fPvslambda.close();
            }
        }

// ************************************************************************

void Output_PvsE(bool lasttime,ofstream &fPvsE,string outputfilenamestem,double lambda,vector<vector<array<double,NY> > > &Y,vector<vector<MATRIX<complex<double>,NF,NF> > > &Scumulative, bool ecsvformat)
      { string cmdotdat("cm.dat");
        stringstream filename;

        if(lasttime==false){ filename.str(""); filename<<outputfilenamestem<<string(":PvsE:lambda=")<<lambda<<cmdotdat; fPvsE.open((filename.str()).c_str()); fPvsE.precision(12);}
        else{ filename.str(""); filename<<outputfilenamestem<<string(":PvsE:lambda=")<<lambdamax0<<cmdotdat; fPvsE.open((filename.str()).c_str()); fPvsE.precision(12);}        

        double r, rrho, YYe;

        // ******

        r = sqrt( RE*RE + lambda*lambda - 2.*RE*lambda*sin(-altitude) );
        if(r > RE){ r = RE;}              
        rrho=rho(r);
        YYe=Ye(r); 

        array<MATRIX<complex<double>,NF,NF>,NM> VfMSW;

        VfMSW[nu][e][e]=Ve(rrho,YYe);
        VfMSW[nu][mu][mu]=Vmu(rrho,YYe);
        VfMSW[nu][tau][tau]=Vtau(rrho,YYe);
        VfMSW[antinu]=-VfMSW[nu];

        vector<vector<MATRIX<complex<double>,NF,NF> > > Hf(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));
        vector<vector<MATRIX<complex<double>,NF,NF> > > UU(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));        
        vector<vector<MATRIX<complex<double>,NF,NF> > > Sa(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)), Sm(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)), Smf(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)), Sfm(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)), Sf(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)); 

        vector<vector<array<double,NF> > > kk(NM,vector<array<double,NF> >(NE));
        vector<vector<array<double,NF> > > dkk(NM,vector<array<double,NF> >(NE));

        int i;
        #pragma omp parallel for schedule(static)
        for(i=0;i<=NE-1;i++)
           { Hf[nu][i]=HfV[nu][i] + VfMSW[nu];
             kk[nu][i]=k(Hf[nu][i]);
               dkk[nu][i]=deltak(kk[nu][i]);
             UU[nu][i]=MixingMatrix(Hf[nu][i],kk[nu][i],dkk[nu][i]);

             Sa[nu][i] = W(Y[nu][i]) * B(Y[nu][i]);

             if(lasttime==false){ Sm[nu][i] = Sa[nu][i] * Scumulative[nu][i];}
             else{ Sm[nu][i] = Adjoint(UV[nu])*UU[nu][i] * Sa[nu][i] * Scumulative[nu][i];}
             
             Smf[nu][i]= Sm[nu][i] * Adjoint(UV[nu]);
             
             if(lasttime==false){ 
                Sfm[nu][i] = UU[nu][i] * Sm[nu][i];
                Sf[nu][i] = UU[nu][i] * Smf[nu][i];
               } 
             else{ Sfm[nu][i] = UV[nu] * Sm[nu][i];
                   Sf[nu][i] = UV[nu] * Smf[nu][i];
                  }

             // *********
             
             Hf[antinu][i]=HfV[antinu][i] + VfMSW[antinu];
               kk[antinu][i]=kbar(Hf[antinu][i]);
             dkk[antinu][i]=deltakbar(kk[antinu][i]);
             UU[antinu][i]=MixingMatrix(Hf[antinu][i],kk[antinu][i],dkk[antinu][i]);
       
             Sa[antinu][i] = W(Y[antinu][i]) * B(Y[antinu][i]);

             if(lasttime==false){ Sm[antinu][i] = Sa[antinu][i] * Scumulative[antinu][i];}
             else{ Sm[antinu][i] = Adjoint(UV[antinu])*UU[antinu][i] * Sa[antinu][i] * Scumulative[antinu][i];}
                          
             Smf[antinu][i]= Sm[antinu][i] * Adjoint(UV[antinu]);
             
             if(lasttime==false)
               { Sfm[antinu][i] = UU[antinu][i] * Sm[antinu][i];
                     Sf[antinu][i] = UU[antinu][i] * Smf[antinu][i];
                    }
             else{ Sfm[antinu][i] = UV[antinu] * Sm[antinu][i];
                       Sf[antinu][i] = UV[antinu] * Smf[antinu][i];
                      }
            }

        // *******
        
        // *******

        fPvsE<<"E [MeV]";

        fPvsE<<"\t P11 \t P12 \t P13 \t P21 \t P22 \t P23 \t P31 \t P32 \t P33";
        fPvsE<<"\t Pbar11 \t Pbar12 \t Pbar13 \t Pbar21 \t Pbar22 \t Pbar23 \t Pbar31 \t Pbar32 \t Pbar33";

        fPvsE<<"\t Pe1 \t Pe2 \t Pe3 \t Pmu1 \t Pmu2 \t Pmu3 \t Ptau1 \t Ptau2 \t Ptau3";
        fPvsE<<"\t Pbare1 \t Pbare2 \t Pbare3 \t Pbarmu1 \t Pbarmu2 \t Pbarmu3 \t Pbartau1 \t Pbartau2 \t Pbartau3";        

        fPvsE<<"\t Pee \t Pemu \t Petau \t Pmue \t Pmumu \t Pmutau \t Ptaue \t Ptaumu \t Ptautau";
        fPvsE<<"\t Pbaree \t Pbaremu \t Pbaretau \t Pbarmue \t Pbarmumu \t Pbarmutau \t Pbartaue \t Pbartaumu \t Pbartautau";        

        for(i=0;i<=NE-1;i++)
           { fPvsE<<"\n"<<E[i]/(mega*cgs::units::eV); 

             fPvsE<<"\t"<<norm(Sm[nu][i][0][0])<<"\t"<<norm(Sm[nu][i][0][1])<<"\t"<<norm(Sm[nu][i][0][2]);
             fPvsE<<"\t"<<norm(Sm[nu][i][1][0])<<"\t"<<norm(Sm[nu][i][1][1])<<"\t"<<norm(Sm[nu][i][1][2]);
             fPvsE<<"\t"<<norm(Sm[nu][i][2][0])<<"\t"<<norm(Sm[nu][i][2][1])<<"\t"<<norm(Sm[nu][i][2][2]);

             fPvsE<<"\t"<<norm(Sm[antinu][i][0][0])<<"\t"<<norm(Sm[antinu][i][0][1])<<"\t"<<norm(Sm[antinu][i][0][2]);
             fPvsE<<"\t"<<norm(Sm[antinu][i][1][0])<<"\t"<<norm(Sm[antinu][i][1][1])<<"\t"<<norm(Sm[antinu][i][1][2]);
             fPvsE<<"\t"<<norm(Sm[antinu][i][2][0])<<"\t"<<norm(Sm[antinu][i][2][1])<<"\t"<<norm(Sm[antinu][i][2][2]);

             fPvsE<<"\t"<<norm(Sfm[nu][i][e][0])<<"\t"<<norm(Sfm[nu][i][e][1])<<"\t"<<norm(Sfm[nu][i][e][2]);
             fPvsE<<"\t"<<norm(Sfm[nu][i][mu][0])<<"\t"<<norm(Sfm[nu][i][mu][1])<<"\t"<<norm(Sfm[nu][i][mu][2]);
             fPvsE<<"\t"<<norm(Sfm[nu][i][tau][0])<<"\t"<<norm(Sfm[nu][i][tau][1])<<"\t"<<norm(Sfm[nu][i][tau][2]);

             fPvsE<<"\t"<<norm(Sfm[antinu][i][e][0])<<"\t"<<norm(Sfm[antinu][i][e][1])<<"\t"<<norm(Sfm[antinu][i][e][2]);
             fPvsE<<"\t"<<norm(Sfm[antinu][i][mu][0])<<"\t"<<norm(Sfm[antinu][i][mu][1])<<"\t"<<norm(Sfm[antinu][i][mu][2]);
             fPvsE<<"\t"<<norm(Sfm[antinu][i][tau][0])<<"\t"<<norm(Sfm[antinu][i][tau][1])<<"\t"<<norm(Sfm[antinu][i][tau][2]);

             fPvsE<<"\t"<<norm(Sf[nu][i][e][e])<<"\t"<<norm(Sf[nu][i][e][mu])<<"\t"<<norm(Sf[nu][i][e][tau]);
             fPvsE<<"\t"<<norm(Sf[nu][i][mu][e])<<"\t"<<norm(Sf[nu][i][mu][mu])<<"\t"<<norm(Sf[nu][i][mu][tau]);
             fPvsE<<"\t"<<norm(Sf[nu][i][tau][e])<<"\t"<<norm(Sf[nu][i][tau][mu])<<"\t"<<norm(Sf[nu][i][tau][tau]);

             fPvsE<<"\t"<<norm(Sf[antinu][i][e][e])<<"\t"<<norm(Sf[antinu][i][e][mu])<<"\t"<<norm(Sf[antinu][i][e][tau]);
             fPvsE<<"\t"<<norm(Sf[antinu][i][mu][e])<<"\t"<<norm(Sf[antinu][i][mu][mu])<<"\t"<<norm(Sf[antinu][i][mu][tau]);
             fPvsE<<"\t"<<norm(Sf[antinu][i][tau][e])<<"\t"<<norm(Sf[antinu][i][tau][mu])<<"\t"<<norm(Sf[antinu][i][tau][tau]);
            }

         fPvsE.flush();
         fPvsE.close();
        }

// ************************************************************************

void Output_Hvslambda(bool firsttime,bool lasttime,ofstream &fHvslambda,double lambda,vector<vector<array<double,NY> > > &Y,vector<vector<MATRIX<complex<double>,NF,NF> > > &Scumulative, bool ecsvformat)
          { MATRIX<complex<double>,NF,NF> VfMSW,VfMSWbar;
            double r, rrho, YYe;

            // *************

            r = sqrt( RE*RE + lambda*lambda - 2.*RE*lambda*sin(-altitude) );
            if(r > RE){ r = RE;}                  
            rrho=YYe=0.; 

            // ****************

            if(firsttime==false && lasttime==false)
              { rrho=rho(r);
                YYe=Ye(r); 
            
                VfMSW[e][e]=Ve(rrho,YYe); 
                VfMSW[mu][mu]=Vmu(rrho,YYe);
                VfMSW[tau][tau]=Vtau(rrho,YYe);
                
                VfMSWbar=-Conjugate(VfMSW);                
               } 

            // **************
         
            if(firsttime==true){
               fHvslambda<<"lambda [cm] \t r [cm] \t rho [g/cm^3] \t Ye [] \t HMSW_ee [erg]";
              }         

            if(firsttime==true){ fHvslambda<<"\n"<<lambdamin0;}
            if(lasttime==true){ fHvslambda<<"\n"<<lambdamax0;}            
            if(firsttime==false && lasttime==false){ fHvslambda<<"\n"<<lambda;}            
            
            fHvslambda<<"\t"<<r<<"\t"<<rrho<<"\t"<<YYe<<"\t"<<real(VfMSW[e][e]);

            fHvslambda.flush();
           }

