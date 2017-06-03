// $Id: jack_corrfunc_ratio.cc,v 2.0 2008/12/05 04:44:06 edwards Exp $
// $Log: jack_corrfunc_ratio.cc,v $
// Revision 2.0  2008/12/05 04:44:06  edwards
// Changed to version 2.0.
//
// Revision 1.2  2007/04/19 01:50:12  kostas
// fixed bug
//
// Revision 1.1  2007/04/18 21:22:31  kostas
// added an other beauty
//
//

#include <iostream>
#include <cstdio>

#include <covfit/statistics.h>
#include <covfit/Function.h>
#include <covfit/fitters/fitfunctions.h>
#include <covfit/fitter.h>
#include <covfit/jackknife.h>

//#include "jackknife_io.h"
//#include "basic_functions.h"

using namespace std;
using namespace CovFit ;


void write(const std::string& t,
           const Array<double>& v,
           const Array<double>& e){
  //Needs size checking                                                        
  for(int i(0);i<v.size();i++)
    cout<<t<<": "<<i<<" "<<v[i]<<" "<<e[i]<<endl ;
}

int main(int argc, char **argv)
{
  if ( (argc != 4 ) && ( argc !=5 ) )
  {
    cerr << "Usage: " << argv[0] << " ";
    cerr<< "<ratio proplist> <two particle proplist> <particle 1 proplist> ";
    cerr<< "[<particle 2 proplist>]" ;
    cerr << endl;
    exit(1);
  }

  PropList twoP ;
  PropList P1,P2 ;

  int Nx = ReadProplist(twoP,argv[2]);
  ReadProplist(P1,argv[3]);
  // if P2 is not present use P1
  if(argc == 5 )
    ReadProplist(P2,argv[4]);
  else
    P2 = P1 ;

  

  
  int Ncnfs(twoP.size()) ;
  int Nt(twoP[0].size());
  if ( ( Ncnfs != P1.size()) || ( Ncnfs != P2.size())  )
    {
      cerr << argv[0] << " : Incompatible configuration numbers" << endl;
      exit(2);
    }
  if( ( Nt != P1[0].size() ) || ( Nt != P2[0].size() ) )
    {
      cerr << argv[0] << " : Time lengths not equal" << endl;
      exit(3);
    }

  PropList jmres(Ncnfs) ;
  Array<double> fmres(Ncnfs) ;

  //Jackknife loop 
  for(int j(0);j<Ncnfs;j++)
    {
      PropList jtwoP ;
      PropList jP1,jP2 ;
      splice(jtwoP,twoP,j,j) ;
      splice(jP1,P1,j,j) ;
      splice(jP2,P2,j,j) ;
      Array<double> mtwoP = mean(jtwoP) ;
      Array<double> mpp = mean(jP1)*mean(jP2) ;
      jmres[j] = mtwoP/mpp ;
    }

  Array<double> mres, e_mres ;
  
  mres = mean(jmres) ;
  e_mres = jackerr(jmres) ;

  write("RATIO",mres,e_mres) ;



  Array<double> mass, e_mass ;
  Array<Array<double> > jmass ;

  jmass = effmass(jmres) ;

  mass = mean(jmass) ;
  e_mass = jackerr(jmass) ;
  
  write("EFFMASS", mass, e_mass) ;

  PropList ratio(Ncnfs);
   //rescale up fluctuations
  for(int j(0);j<Ncnfs;j++)
    {
      Array<double> d = jmres[j] - mres ;
      ratio[j] = mres + double(Ncnfs-1)*d ;
    }
  WriteProplist(ratio,Nx,argv[1]);

}
