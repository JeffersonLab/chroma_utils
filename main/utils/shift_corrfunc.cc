// $Id: shift_corrfunc.cc,v 2.0 2008/12/05 04:44:06 edwards Exp $
// $Log: shift_corrfunc.cc,v $
// Revision 2.0  2008/12/05 04:44:06  edwards
// Changed to version 2.0.
//
// Revision 1.1  2005/10/30 19:15:43  kostas
// Added the correlation function shift utility
//
//

#include <iostream>
#include <cstdio>

#include <covfit/proplist.h>


using namespace std;
using namespace CovFit;

int main(int argc, char **argv)
{
  if ((argc < 4))
  {
    cerr << "Usage: " << argv[0] << " <tsrc> <out proplist> <in proplist>" ; 
    cerr << endl;
    exit(1);
  }

  int tsrc(atoi(argv[1]));
  
  PropList pp ;
  PropList sh ;

  int Nx = ReadProplist(pp,argv[3]);
  int Nt(pp[0].size()); //assume all props have same length 
  cout<< "Time length: "<<Nt<<endl ;
  cout<< "Time shift : "<<tsrc<<endl ;

  sh.resize(pp.size());
  for(int j(0);j<pp.size();j++){
    sh[j].resize(Nt);
    for(int i(0);i<Nt;i++){
      int ii((i-tsrc+Nt)%Nt) ;
      //cout<<"i: "<<i<<" ii: "<<ii<<endl ;
      sh[j][ii] = pp[j][i] ;
    }
  }
    
  WriteProplist(sh,Nx,argv[2]);
}
