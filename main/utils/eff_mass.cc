// eff_mass writing out  a jackknife file
//
//

#include <iostream>
#include <cstdio>

#include "covfit/proplist.h"
#include <covfit/jackknife.h>

using namespace std;
using namespace CovFit ;

void write(const std::string& t,
           const Array<Double>& v,
           const Array<Double>& e){
  //Needs size checking                                                        
  for(int i(0);i<v.size();i++)
    cout<<t<<": "<<i<<" "<<v[i]<<" "<<e[i]<<endl ;
}

void error(int argc, char **argv){
  cerr << "Usage: " << argv[0] << " [-t exp (default)|cosh|pcosh] [-d 1(default) <time split>] <proplist>" ;
  cerr <<endl;
  cerr << "   exp  : exponential (solve two neighboring points)"<< endl;
  cerr << "   cosh : exponential (fit three points)"<< endl;
  cerr << "   stag : exponential (fit three points with step 2)"<< endl;
  cerr << "   [cosh] and [stag] work for periodic BC as well as Dirichlet BC" ;
  cerr << endl;
  cerr << endl<<endl;
  cerr<< "Found "<<argc<<endl ;
  cerr<< "Command line: ";
  for(int i(0);i<argc; i++)
    cerr<<argv[i]<<" ";
  cerr<<endl ;
  exit(1);
}


int main(int argc, char **argv)
{
  if ((argc!=4)&& (argc!=2)&&(argc!=6) ) error(argc,argv) ;

  int type(0);
  int arg_ofset = 0 ;
  if(argc==4) 
    arg_ofset = 2 ;
  if(argc==6) 
    arg_ofset = 4 ;

  int t_arg(0);
  int d_arg(0);

  if(string(argv[1]) == "-t"){
    t_arg = 2 ;
    if(argc < 4 ) error(argc,argv) ;
  }
  if(string(argv[3]) == "-t"){
    t_arg = 4 ;
    if(argc < 6 ) error(argc,argv) ;
  }
  if(string(argv[1]) == "-d"){
    d_arg = 2 ;
    if(argc < 4 ) error(argc,argv) ;
  }
  if(string(argv[3]) == "-d"){
    d_arg = 4 ;
    if(argc < 6 ) error(argc,argv) ;
  }

  if(t_arg){
    //cout<<argv[1]<<endl ;
    if (string(argv[t_arg]) == "cosh" )
      type = 1 ;
    else  if (string(argv[t_arg]) == "stag" )
      type = 2 ;
    else if (string(argv[t_arg]) == "exp" )
      type = 0 ;
    else{
      cerr<<" Unknown option: -t "<<argv[t_arg]<<endl ;
      error(argc,argv) ;
    }
  }
  int d(1) ;
  if(d_arg){
    //cout<<argv[1]<<endl ;
    d=atoi(argv[d_arg]) ;
  }


  PropList pp ;
  PropList jmass ;
  int Nx = ReadProplist(pp,argv[1+arg_ofset]);
  Array< Array <Double> > jpp(pp.size()) ;
  
  //cout<<Nx<<endl ;
  for(int j(0);j<pp.size();j++){
    PropList tt ;
    splice(tt,pp,j,j) ;
    jpp[j] = mean(tt) ;
  }

  cout<<"Using effective mass type: "<< type<<endl ;
  cout<<"Using distance: "<< d<<" to calculate the effective mass"<<endl ;

  switch ( type ) {
  case 1: 
    jmass = effmass_cosh(jpp,d) ;
    break ;
  case 2: 
    if(d!=1) 
      cout<<"WARNING: Staggered effective mass works with d=1 only."<<endl ;
    jmass = effmass_cosh_stag(jpp) ;
    break ;
  case 0:
  default:
    jmass = effmass(jpp,d) ;
  }


  Array<Double> mass = mean(jmass) ;
  Array<Double> e_mass = jackerr(jmass) ;

  write("EFFMASS", mass, e_mass) ;
  

  WriteProplist(jmass,Nx,"jk_eff_mass.lst");

}
