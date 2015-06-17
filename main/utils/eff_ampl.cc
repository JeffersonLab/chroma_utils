// eff_ampl writing out  a jackknife file
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
  cerr << "Usage: " << argv[0] << " [-d 1(default) <time split>] <proplist>" ;
  cerr <<endl;
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
  if ((argc!=4)&& (argc!=2) ) error(argc,argv) ;

  //  int type(0);
  int arg_ofset = 0 ;
  if(argc==4) 
    arg_ofset = 2 ;
  //  if(argc==6) 
  //  arg_ofset = 4 ;

  //int t_arg(0);
  int d_arg(0);

  if(string(argv[1]) == "-d"){
    d_arg = 2 ;
    if(argc < 4 ) error(argc,argv) ;
  }

  int d(1) ;
  if(d_arg){
    //cout<<argv[1]<<endl ;
    d=atoi(argv[d_arg]) ;
  }


  PropList pp ;
  int Nx = ReadProplist(pp,argv[1+arg_ofset]);
  Array< Array <Double> > jpp(pp.size()) ;
  PropList jAmp(pp.size()) ;

  cout<<"Nx = "<<Nx<<endl ;
  for(int j(0);j<pp.size();j++){
    PropList tt ;
    splice(tt,pp,j,j) ;
    jpp[j] = mean(tt) ;
  }

  //cout<<"Using effective mass type: "<< type<<endl ;
  cout<<"Using distance: "<< d<<" to calculate the effective amplitude"<<endl ;

  int Nt=pp[0].size();
  cout<<"Nt = "<<Nt<<endl ;
  for(int j(0);j<pp.size();j++){
    cout<<"Doing jacknife sample: "<<j<<endl ;
    Array<Double> jj(Nt);
    for(int t(0);t<Nt;t++){
      if((t+d)<Nt)
	jj[t] = jpp[j][t]*pow(jpp[j][t]/jpp[j][t+d],double(t)/double(d));
      else
	jj[t] = jj[Nt-d-1] ;
      //cout<<" jj["<<t<<"]="<<jj[t]<<endl ;
    }
    jAmp[j]=jj ;
    //cout<<"jAmp.size()="<<jAmp.size()<<endl ;
  }


  Array<Double> Ampl = mean(jAmp) ;
  Array<Double> e_Ampl = jackerr(jAmp) ;

  write("EFFAMPL", Ampl, e_Ampl) ;
  

  WriteProplist(jAmp,Nx,"jk_eff_ampl.lst");

}
