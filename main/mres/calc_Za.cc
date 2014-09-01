// $Id: calc_Za.cc,v 2.0 2008/12/05 04:43:50 edwards Exp $
// $Log: calc_Za.cc,v $
// Revision 2.0  2008/12/05 04:43:50  edwards
// Changed to version 2.0.
//
// Revision 1.5  2005/11/23 04:32:10  kostas
// ..
//
// Revision 1.4  2005/04/27 19:53:52  edwards
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
// Revision 1.3  2004/11/19 23:59:00  kostas
// fixed calc_Za
//
// Revision 1.2  2004/03/17 16:29:34  kostas
// Za calculation; dwf locality stripper
//
// Revision 1.1  2004/03/12 16:29:42  kostas
// Za calculation
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
    cerr << "Usage: " << argv[0] << " <conAx proplist> <locAx proplist> [<tmin>, <tmax>]" ; 
    cerr << endl;
    exit(1);
  }

  PropList conAx ;
  PropList locAx ;

  ReadProplist(conAx,argv[1]);
  ReadProplist(locAx,argv[2]);

  
  int Ncnfs(conAx.size()) ;
  int Nt(conAx[0].size());
  if( Ncnfs != locAx.size())
    {
      cerr << argv[0] << " : Incompatible configuration numbers" << endl;
      exit(2);
    }
  if( Nt != locAx[0].size() )
    {
      cerr << argv[0] << " : Time lengths not equal" << endl;
      exit(3);
    }

  PropList jZa(Ncnfs) ;
  for(int j(0);j<Ncnfs;j++)
    jZa[j].resize(Nt) ;

  Array<Double> fZa(Ncnfs) ;
  Polynomium P(1);
  FitterArgs fargs ;
  
  fargs.xlow = (4) ;
  fargs.xhigh = (17);
  if(argc >= 5 ){
    fargs.xlow  = atoi(argv[3]);
    fargs.xhigh = atoi(argv[4]);
  }
  fargs.npar = P.Npar();
  Fitter fit(P,fargs);
  Array<Double> time(Nt) ;
  for(int t(0);t<Nt;t++)
    time[t] = t ;

  //Jackknife loop 
  for(int j(0);j<Ncnfs;j++)
    {
      PropList jconAx ;
      PropList jlocAx ;
      splice(jconAx,conAx,j,j) ;
      splice(jlocAx,locAx,j,j) ;
      Array<Double> mconAx = mean(jconAx) ;
      Array<Double> econAx = err (jconAx) ;
      Array<Double> mlocAx = mean(jlocAx) ;
      Array<Double> elocAx = err (jlocAx) ;
      Array<Double> err(Nt) ;
      jZa[j][0] = 2.0*mconAx[0]/(mlocAx[0] + mlocAx[1]) ;
      jZa[j][Nt-1] = (mconAx[Nt-1] + mconAx[Nt-2])/(2.0*mlocAx[Nt-1]) ;
      // not correct but who cares for point 0 and Nt-1
      err[0] = econAx[0] + elocAx[0] + elocAx[1] ; 
      err[Nt-1] =econAx[Nt-1] + econAx[Nt-2] + elocAx[Nt-1]   ;

      for(int t(1);t<Nt-1;t++){
	Double t1( (mconAx[t] + mconAx[t-1])/(2.0*mlocAx[t]) );
	Double t2( 2.0*mconAx[t]/(mlocAx[t] + mlocAx[t+1])   );
	jZa[j][t] = 0.5*( t1 + t2 ) ;
	Double elt(elocAx[t]*(t1/mlocAx[t]+2.0*mconAx[t]/mlocAx[t]/mlocAx[t]));
	Double eltt(elocAx[t+1]*(2.0*mconAx[t]/mlocAx[t+1]/mlocAx[t+1]));
	Double ect(econAx[t]*(0.5/mlocAx[t]+t2/mconAx[t]));
	Double ectt(econAx[t+1]*0.5/mlocAx[t]);
       
	err[t] = sqrt(elt*elt + eltt*eltt + ect*ect+ ectt*ectt);
      }

      Array<Double> v(1),e(1);
      if(j==0)
	v[0] = .8 ;
      else
	v[0] = fZa[j-1] ;
      fit.setParams(v) ;
      fit.setData(time,jZa[j],err) ;
      fit.Fit();
      fit.getParams(v, e) ;
      fZa[j] = v[0] ;
    }
  
  Array<Double> Za, e_Za ;
  
  Za = mean(jZa) ;
  e_Za = jackerr(jZa) ;

  for(int t(0);t<Nt;t++)
    cout<<"ZA: "<<t<<" "<<Za[t]<<" "<< e_Za[t]<<endl ;

  Double val, e_val ;
  val = mean(fZa) ;
  e_val = jackerr(fZa) ;

  cout<<"Jackknife fitted Za: "<< val << " +/- " <<e_val<<endl ;
}
