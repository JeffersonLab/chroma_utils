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
  if (argc != 4)
  {
    cerr << "Usage: " << argv[0] << " tmax <out proplist> <in proplist>" ; 
    cerr << endl;
    cerr << "   prints proplist from t=0 to t=tmax-1 \n" ;
    exit(1);
  }

  PropList p ;
  

  int Nx=ReadProplist(p,argv[3]);

  
  int Ncnfs(p.size()) ;
  int Nt(p[0].size());
  int tmax = atoi(argv[1]) ;
  if((tmax>Nt)||(tmax<0)){
    cerr<<" tmax is out of range"<<endl ;
    exit(2);
  }
  PropList cp(Ncnfs) ;
  for(int j(0);j<Ncnfs;j++){
    cp[j].resize(tmax) ;
    for(int t(0);t<tmax;t++)
      cp[j][t] = p[j][t];
  }

  WriteProplist(cp,Nx,argv[2]);
    
}
