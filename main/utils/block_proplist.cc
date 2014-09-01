// $Id: block_proplist.cc,v 2.0 2008/12/05 04:44:06 edwards Exp $
// $Log: block_proplist.cc,v $
// Revision 2.0  2008/12/05 04:44:06  edwards
// Changed to version 2.0.
//
// Revision 1.2  2007/10/22 03:30:21  kostas
// fixed a bug in block utility
//
// Revision 1.1  2006/10/05 14:04:33  kostas
// added blocking utility program
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
#include "covfit/proplist.h"

using namespace std;
using namespace CovFit;

int main(int argc, char **argv)
{
  if (argc != 4)
  {
    cerr << "Usage: " << argv[0] << " <proplist> <independent cnf number> <out proplist>" ; 
    cerr << endl;
    exit(1);
  }

  PropList in ;

  int Nx = ReadProplist(in,argv[1]);
  
  int cnfs ;
  cnfs = atoi(argv[2]);

  
  int N(in.size()) ;
  int block(N/cnfs);
  cout<<"Block factor is "<<Real(N)/Real(cnfs)<<endl;
  int Nt(in[0].size());


  PropList out(cnfs) ;


  //Jackknife loop 
  for(int j(0);j<cnfs;j++)
    {
      PropList blk ;
      arrayblock(blk,in,j*block,(j+1)*block-1) ;
      out[j] = mean(blk) ;
    }
  
  WriteProplist( out, Nx, argv[3]);
}
