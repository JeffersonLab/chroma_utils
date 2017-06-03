// $Id: jack_ratio.cc,v 2.0 2008/12/05 04:44:06 edwards Exp $
// $Log: jack_ratio.cc,v $
// Revision 2.0  2008/12/05 04:44:06  edwards
// Changed to version 2.0.
//
// Revision 1.2  2005/11/23 03:16:22  kostas
// moved fitter headers into the main include tree
//
// Revision 1.1  2005/10/30 21:17:25  kostas
// added jack_ratio utility
//
//
//

#include <iostream>
#include <cstdio>

#include <covfit/statistics.h>
#include <covfit/Function.h>
#include <covfit/fitters/fitfunctions.h>
#include <covfit/fitter.h>

using namespace std;
using namespace CovFit;

int main(int argc, char **argv)
{
  if (argc != 5)
  {
    cerr << "Usage: " << argv[0] << " <fit tmin> <fit tmax> <numerator proplist> <denumerator proplist>" ; 
    cerr << endl;
    cerr << "   jackknifes a ratio and  fits a constant between tmax and tmin\n" ;
    exit(1);
  }

  PropList mp ;
  PropList pp ;

  ReadProplist(mp,argv[3]);
  ReadProplist(pp,argv[4]);

  
  int Ncnfs(mp.size()) ;
  int Nt(mp[0].size());
  if( Ncnfs != pp.size())
    {
      cerr << argv[0] << " : Incompatible configuration numbers" << endl;
      exit(2);
    }
  if( Nt != pp[0].size() )
    {
      cerr << argv[0] << " : Time lengths not equal" << endl;
      exit(3);
    }

  PropList jmres(Ncnfs) ;
  Array<double> fmres(Ncnfs) ;
  Polynomium P(1);
  
  FitterArgs fargs ;
  
  fargs.xlow = atoi(argv[1]) ;
  fargs.xhigh = atoi(argv[2]);
  fargs.npar = P.Npar();
  Fitter fit(P,fargs);
  Array<double> time(Nt) ;
  for(int t(0);t<Nt;t++)
    time[t] = t ;

  Array<double> v(1),e(1);
  //Jackknife loop 
  for(int j(0);j<Ncnfs;j++)
    {
      PropList jmp ;
      PropList jpp ;
      splice(jmp,mp,j,j) ;
      splice(jpp,pp,j,j) ;
      Array<double> mmp = mean(jmp) ;
      Array<double> emp = err (jmp) ;
      Array<double> mpp = mean(jpp) ;
      Array<double> epp = err (jpp) ;
      jmres[j] = mmp/mpp ;
      Array<double> err = jmres[j]*sqrt(emp*emp/(mmp*mmp) + epp*epp/(mpp*mpp));
     
      if(j==0){
	v[0] = jmres[j][1] ;
      }
   
      fit.setParams(v) ;
      fit.setData(time,jmres[j],err) ;
      fit.Fit();
      fit.getParams(v, e) ;
      fmres[j] = v[0] ;
    }
  
  Array<double> mres, e_mres ;
  
  mres = mean(jmres) ;
  e_mres = jackerr(jmres) ;

  for(int t(0);t<Nt;t++)
    cout<<"RATIO: "<<t<<" "<<mres[t]<<" "<< e_mres[t]<<endl ;

  double val, e_val ;
  val = mean(fmres) ;
  e_val = jackerr(fmres) ;

  cout<<"Jackknife fitted RATIO: "<< val << " +/- " <<e_val<<endl ;
}
