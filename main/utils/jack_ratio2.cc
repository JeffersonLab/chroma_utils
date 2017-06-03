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
  if (argc != 3)
  {
    cerr << "Usage: " << argv[0] << " <numerator proplist> <denumerator proplist>" ; 
    cerr << endl;
    cerr << "   jackknifes the squared ratio. Assumes jackknife lists\n" ;
    cerr << "   Writes out both rescaled and regular jackknife list\n" ;
    exit(1);
  }

  PropList mp ;
  PropList pp ;

  int Nx = ReadProplist(mp,argv[1]);
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

  //Jackknife loop 
  // already jackknife lists
  for(int j(0);j<Ncnfs;j++)
    {
      jmres[j] = mp[j]/pp[j] ;
      jmres[j] *= jmres[j] ;
    }
  

  Array<double> mres,e_mres;

  mres = mean(jmres) ;
  e_mres = jackerr(jmres) ;

  for(int t(0);t<Nt;t++)
    cout<<"RATIO2: "<<t<<" "<<mres[t]<<" "<< e_mres[t]<<endl ;

  WriteProplist(jmres,Nx,"jk_ratio2.lst");

  for(int j(0);j<Ncnfs;j++){
    Array<double> d = jmres[j] - mres ;
    jmres[j] = mres + double(Ncnfs-1)*d ;
  }
  WriteProplist(jmres,Nx,"ratio2.lst");
}
