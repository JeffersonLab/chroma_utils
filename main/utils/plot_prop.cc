// $Id: plot_prop.cc,v 2.0 2008/12/05 04:44:06 edwards Exp $
// $Log: plot_prop.cc,v $
// Revision 2.0  2008/12/05 04:44:06  edwards
// Changed to version 2.0.
//
// Revision 1.1  2007/05/10 19:58:25  kostas
// added some utils
//
// Revision 1.3  2005/12/30 02:21:22  kostas
// fixed more command line parsing bugs...
//
// Revision 1.2  2005/12/29 17:30:40  kostas
// fixed antiperiodic command line bug
//
// Revision 1.1  2005/10/14 01:33:36  kostas
// added propagator fold code
//
//

#include <iostream>
#include <cstdio>

#include "covfit/proplist.h"
#include "covfit/statistics.h"

using namespace std;
using namespace CovFit ;

void error(int argc, char **argv){
  cerr << "Usage: " << argv[0] << "  <proplist>\n";
  cerr<< "Found "<<argc<<endl ;
  cerr<< "Command line: ";
  for(int i(0);i<argc; i++)
    cerr<<argv[i]<<" ";
  cerr<<endl ;
  exit(1);
}


int main(int argc, char **argv)
{
  if (argc != 2 ) error(argc,argv) ;

  

  PropList pp ;
  Array<double> av,er ;
  int Nx = ReadProplist(pp,argv[1]);

  av = mean(pp) ;
  //err = jackerr(pp) ;
  er = err(pp) ;

  for(int t(0);t<av.size();t++)
    cout<<t<<" "<<av[t]<<" "<<er[t]<<endl ;
  
}
