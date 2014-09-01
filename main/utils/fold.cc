// $Id: fold.cc,v 2.0 2008/12/05 04:44:06 edwards Exp $
// $Log: fold.cc,v $
// Revision 2.0  2008/12/05 04:44:06  edwards
// Changed to version 2.0.
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

using namespace std;
using namespace CovFit ;

void error(int argc, char **argv){
  cerr << "Usage: " << argv[0] << " [<-bc b>] <TIME REVERSED proplist> <proplist>";
  cerr << "\n     b is 1 for periodic (default) -1 for antiperiodic"<< endl;
  cerr<< "Found "<<argc<<endl ;
  cerr<< "Command line: ";
  for(int i(0);i<argc; i++)
    cerr<<argv[i]<<" ";
  cerr<<endl ;
  exit(1);
}


int main(int argc, char **argv)
{
  if ((argc != 3) && ( argc != 5) ) error(argc,argv) ;

  int arg_ofset = 0 ;
  double sign(1.0);
  if(string(argv[1]) == "-bc"){
    arg_ofset = 2 ;
    sign = atof(argv[2]) ;
    if(argc != 5 ) error(argc,argv) ;
  }


  PropList pp ;
  PropList av ;
  int Nx = ReadProplist(pp,argv[2+arg_ofset]);

  av.resize(pp.size());
  for(int j(0);j<pp.size();j++){
    av[j].resize(pp[j].size()/2+1); // assume even pp[j].size() 
    av[j][0] = pp[j][0] ;
    for(int i(1);i<av[j].size();i++){
      av[j][i]  = .5*(pp[j][i] + sign*pp[j][pp[j].size()-i]) ;
    }
  }
	     
  WriteProplist(av,Nx,argv[1+arg_ofset]);
  
}
