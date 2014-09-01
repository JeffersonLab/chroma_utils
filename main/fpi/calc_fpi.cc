// $Id: calc_fpi.cc,v 2.0 2008/12/05 04:43:50 edwards Exp $
// $Log: calc_fpi.cc,v $
// Revision 2.0  2008/12/05 04:43:50  edwards
// Changed to version 2.0.
//
// Revision 1.5  2005/11/23 04:36:59  kostas
//
// ..
//
// Revision 1.4  2005/11/22 03:00:11  kostas
// adjusted fpi so that it can do fK
//
// Revision 1.3  2005/07/29 14:36:17  edwards
// Changed includes so old adat files are now down in "io" subdir.
//
// Revision 1.2  2005/05/13 19:27:50  kostas
// added fpi main program
//
// Revision 1.1  2005/05/13 16:55:19  kostas
// added fpi
//

#include <iostream>
#include <cstdio>

#include "covfit/statistics.h"
#include "covfit/Function.h"
#include "covfit/fitters/fitfunctions.h"
#include "covfit/fitter.h"
#include "covfit/FpiParams_io.h"

using namespace std;
using namespace CovFit;

int main(int argc, char **argv)
{

  if (argc != 2)
  {
    cerr << "Usage: " << argv[0] << " <params>" ; 
    cerr << endl;
    exit(1);
  }

  PropList corrSP ;
  PropList corrSS ;

  FpiParams_t param ;
  
  
  XMLReader xml_in(argv[1]);
  try{
    read(xml_in,"/FpiParam_db", param) ;
  }
  catch(const string& error) {
    cerr<<"OOOPS!! : "<<error<<"\n" ;
    exit(22);
  }
  xml_in.close();
  
  // ofstream ftl("fpi-jack-list");


  ReadProplist(corrSP,param.filenameSP);
  ReadProplist(corrSS,param.filenameSS);
    
  int Ncnfs(corrSP.size()) ;
  int tt(corrSS.size()) ;
  if(tt!=Ncnfs){
     cerr<<"OOOPS!! : Configuration number missmatch!\n" ;
     exit(1);
  }

  int Nt(corrSP[0].size());
  tt =(corrSS[0].size()) ;
  if(tt!=Nt){
     cerr<<"OOOPS!! : Time legth missmatch!\n" ;
     exit(2);
  }
  
  // For Dirichlet BC simple exponetial fits
  // needs improvement so that different BC can be used
  Exponentials expo(1);
  FitterArgs e_args ;
    
  
  e_args.maxit = param.fit.Niter ;
  e_args.eps   = param.fit.Toler ;
  e_args.npar = expo.Npar();
  e_args.fixpar.resize(0) ; 

  Array<Double> time(Nt) ;
  for(int t(0);t<Nt;t++) time[t] = t ;
    
  int Jblocks = Ncnfs/param.block ;

  int ib,ie ;
  Array<Double> jFpi(Jblocks) ;
  cout<<"Number of configurations: "<<Ncnfs <<endl ;
  cout<<"Number of Jackknife blocks: "<<Jblocks<<endl ;

  for(int j(0);j<Jblocks;j++){
    ib=j*param.block ;
    ie=ib+param.block-1 ;
    PropList jcorrSP, jcorrSS ;

    cout<< "splicing out from "<<ib<<" to "<< ie<<endl ;
    splice(jcorrSP, corrSP, ib,ie);
    splice(jcorrSS, corrSS, ib,ie);


    CovarMat cmSP(Nt,param.cov.block);
    cmSP.SetDataList(time,jcorrSP);
    cmSP.CalcCovarMat();
    cmSP.WriteCovarMat("covarSP.PION");

    CovarMat cmSS(Nt,param.cov.block);
    cmSS.SetDataList(time,jcorrSS);
    cmSS.CalcCovarMat();
    cmSS.WriteCovarMat("covarSS.PION");
    
    
    // CONTINUE with the fitting code
    coFitter efit(expo,e_args);
    Array<Double> guessparSP, fitparSP, fiterrSP ;
    guessparSP = param.fit_paramsSP ;

    Array<Double> guessparSS, fitparSS, fiterrSS ;
    guessparSS = param.fit_paramsSS ;

    //ftl<<"\n"<<param.sta[s].name<<" FITS"<<endl ;
    cout<<"\n\n****************\n" ;  
    cout<<efit.Ndof()<<" PARAMETERS: DISTANCE "<<param.mint<<" to "<<param.maxt<<endl ;

    efit.setRange(param.mint,param.maxt) ;
    efit.setParams(guessparSP);
    efit.SetCovarMat(cmSP) ;
    efit.Fit() ;
    efit.getParams(fitparSP,fiterrSP);

    //efit.setRange(param.mint,param.maxt) ;
    guessparSS[1] = fitparSP[1] ;
    efit.setParams(guessparSS);
    efit.SetCovarMat(cmSS) ;
    efit.Fit() ;
    efit.getParams(fitparSS,fiterrSS);
    
    Double mq (0.5*(param.quark_mass[0]+param.quark_mass[1]) + param.mres) ;
    Double two_over_mpi = 2.0/fitparSP[1] ;
    jFpi(j) = fitparSP[0]/sqrt(fitparSS[0])*mq*two_over_mpi*sqrt(two_over_mpi);

    cout<<"JACK: "<< j << " ";
    cout<<fitparSP[0]<<" +/- "<<fiterrSP[0]<<" " ;
    cout<<fitparSP[1]<<" +/- "<<fiterrSP[1]<<" " ;
    cout<<fitparSS[0]<<" +/- "<<fiterrSS[0]<<" " ;
    cout<<fitparSS[1]<<" +/- "<<fiterrSS[1]<<" " ;
    cout<<jFpi[j] ;
    if(!efit.Convergence()) cout<<" * NO CONVERGENCE: "<<j ;
    if(efit.FailedFit()) cout<<" ! FIT FAILED: "<<j ;
    cout<<endl ;

  }

  Double Fpi = mean(jFpi);
  Double eFpi = jackerr(jFpi) ;
 
  cout<<"Jackknife Fpi: "<<Fpi<<" +/- "<<eFpi<<endl ;
}
