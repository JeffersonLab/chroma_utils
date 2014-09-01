// $Id: mres.cc,v 2.0 2008/12/05 04:43:50 edwards Exp $
// $Log: mres.cc,v $
// Revision 2.0  2008/12/05 04:43:50  edwards
// Changed to version 2.0.
//
// Revision 1.6  2005/11/23 04:32:10  kostas
// ..
//
// Revision 1.5  2005/04/27 19:53:52  edwards
// Completely overhauled include and lib structure. Put all
// of Kostas' fit stuff into their own covfit subdir in include
// and lib. All in its own namespace. Similarly, all of
// ensem is in it's own dirs and namespaces. Unified use of
// XML Array as base array type. I've not touched TTT which
// uses someother classes.
// The new convention should be is people put there little mini
// packages in their own dirs and namespaces. One can compose
// a bigger namespace with namespace composition.
//
// Revision 1.4  2004/11/19 23:59:00  kostas
// fixed calc_Za
//
// Revision 1.3  2004/10/21 21:12:19  kostas
// fixed stripper to read several output versions
//
// Revision 1.2  2004/02/26 20:38:50  kostas
// added jackknifed mres calculation
//
// Revision 1.1  2004/02/10 22:21:06  kostas
// Added some statistics routines an mres computation and some operators
// for the Array class
//

#include <iostream>
#include <cstdio>

#include "covfit/statistics.h"
#include "covfit/Function.h"
#include "covfit/fitters/polynomium.h"
#include "covfit/fitter.h"

using namespace CovFit;
using namespace std;

int main(int argc, char **argv)
{
  if (argc < 3)
  {
    cerr << "Usage: " << argv[0] << " <Midpoint proplist> <Psuedo proplist> [<tmin>, <tmax>]" ; 
    cerr << endl;
    exit(1);
  }

  PropList mp ;
  PropList pp ;

  ReadProplist(mp,argv[1]);
  ReadProplist(pp,argv[2]);

  
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
  Array<Double> fmres(Ncnfs) ;
  Polynomium P(1);
  FitterArgs fargs ;
  
  fargs.xlow = (0) ;
  fargs.xhigh = (17);
  if(argc >= 5 ){
    fargs.xlow  = atoi(argv[3]);
    fargs.xhigh = atoi(argv[4]);
  }

  cout << "Fitting range: "<<fargs.xlow <<" "<<fargs.xhigh<<endl ;
  fargs.npar = P.Npar();
  Fitter fit(P,fargs);
  Array<Double> time(Nt) ;
  for(int t(0);t<Nt;t++)
    time[t] = t ;

  //Jackknife loop 
  for(int j(0);j<Ncnfs;j++)
    {
      PropList jmp ;
      PropList jpp ;
      splice(jmp,mp,j,j) ;
      splice(jpp,pp,j,j) ;
      Array<Double> mmp = mean(jmp) ;
      Array<Double> emp = err (jmp) ;
      Array<Double> mpp = mean(jpp) ;
      Array<Double> epp = err (jpp) ;
      jmres[j] = mmp/mpp ;
      Array<Double> err = jmres[j]*sqrt(emp*emp/(mmp*mmp) + epp*epp/(mpp*mpp));
      Array<Double> v(1),e(1);
      if(j==0)
	v[0] = 1.0e-3 ;
      else
	v[0] = fmres[j-1] ;
      fit.setParams(v) ;
      fit.setData(time,jmres[j],err) ;
      fit.Fit();
      fit.getParams(v, e) ;
      fmres[j] = v[0] ;
    }
  
  Array<Double> mres, e_mres ;
  
  mres = mean(jmres) ;
  e_mres = jackerr(jmres) ;

  for(int t(0);t<Nt;t++)
    cout<<"MRES: "<<t<<" "<<mres[t]<<" "<< e_mres[t]<<endl ;

  Double val, e_val ;
  val = mean(fmres) ;
  e_val = jackerr(fmres) ;

  cout<<"Jackknife fitted residual mass: "<< val << " +/- " <<e_val<<endl ;
}
