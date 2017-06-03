#include "covfit/proplist.h"
#include <iomanip>
#include <fstream>

namespace CovFit
{

  using namespace std;
  
  int ReadProplist(PropList& props,const string& f){
    ifstream file(f.c_str());
    int nprops ;
    int Nt,Nx,fu,bar  ;

    //read the header
    file>>nprops>>Nt>>fu>>Nx>>bar ;
    
    props.resize(nprops);
    for(int i(0);i<props.size();i++)
      props[i].resize(Nt);
    
    double imaginary_part_to_be_discarded ;
    for(int i(0);i<props.size();i++)
      for(int t(0);t<props[i].size();t++)
	{
	  int tt ;
	  file>>tt>>props[i][t] ;
	  if(fu==1) file>>imaginary_part_to_be_discarded ;
	  if(t!=tt){
	    cerr<<"Wrong data in file: "<<f.c_str() ;
	    throw ;
	  }
	}
    
    file.close();

    return Nx ;
  }

  int ReadProplist(cPropList& props,const string& f){
    ifstream file(f.c_str());
    int nprops ;
    int Nt,Nx,fu,bar  ;

    //read the header
    file>>nprops>>Nt>>fu>>Nx>>bar ;
    
    props.resize(nprops);
    for(int i(0);i<props.size();i++)
      props[i].resize(Nt);
    
    for(int i(0);i<props.size();i++)
      for(int t(0);t<props[i].size();t++)
	{
	  int tt ;
	  double re ;
	  file>>tt>>re ; //props[i][t].re ;
	  if(fu==0)
	    props[i][t] = DComplex(re,0.0) ;
	  else{
	    double im ;
	    file>>im ; //props[i][t].im ;
	    props[i][t] = DComplex(re,im) ;
	  }

	  if(t!=tt){
	    cerr<<"Wrong data in file: "<<f.c_str() ;
	    throw ;
	  }
	}
    
    file.close();

    return Nx ;
  }

  void WriteProplist(PropList& props, int Nx, const string& f){
    ofstream file(f.c_str());
    // Write header line
    int nprops(props.size());
    int Nt(props[0].size());
    file<<nprops <<" "<<Nt<<" 0 "<<Nx<<" 1"<<endl;
    file<<setprecision(16) ;
    for(int i(0);i<props.size();i++)
      for(int t(0);t<props[i].size();t++)
	file<<t<<" "<<props[i][t]<<endl ;

    file.close();
  }

  void WriteProplist(cPropList& props, int Nx, const string& f){
    ofstream file(f.c_str());
    // Write header line
    int nprops(props.size());
    int Nt(props[0].size());
    file<<nprops <<" "<<Nt<<" 1 "<<Nx<<" 1"<<endl;
    file<<setprecision(16) ;
    for(int i(0);i<props.size();i++)
      for(int t(0);t<props[i].size();t++){
	file<<t<<" "<<props[i][t].real()<<" "<<props[i][t].imag()<<endl ;
      }
    file.close();
  }

  int ReadProplist(Array<Array< Array<double> > >& props ,const string& f){
    ifstream file(f.c_str());
    int nprops ;
    int Nt,Nx,fu,bar  ;

    //read the header
    file>>nprops>>Nt>>fu>>Nx>>bar ;
    
    props.resize(nprops);
    for(int i(0);i<props.size();i++){
      props[i].resize(Nt);
      for(int k(0);k<props[i].size();k++)
	props[i][k].resize(bar) ;
    }
    
    double imaginary_part_to_be_discarded ;
    for(int i(0);i<props.size();i++)
      for(int t(0);t<props[i].size();t++)
	{
	  int tt ;
	  file>>tt;
	  for(int l(0);l<props[i][t].size();l++){
	    file>>props[i][t][l] ;
	    if(fu==1) file>>imaginary_part_to_be_discarded ;
	  }
	  if(t!=tt){
	    cerr<<"Wrong data in file: "<<f.c_str() ;
	    throw ;
	  }
	}
    
    file.close();

    return Nx ;
  }

  int ReadProplist(Array<Array< Array<DComplex> > >& props ,const string& f){
    ifstream file(f.c_str());
    int nprops ;
    int Nt,Nx,fu,bar  ;

    //read the header
    file>>nprops>>Nt>>fu>>Nx>>bar ;
    
    props.resize(nprops);
    for(int i(0);i<props.size();i++){
      props[i].resize(Nt);
      for(int k(0);k<props[i].size();k++)
	props[i][k].resize(bar) ;
    }
    
    for(int i(0);i<props.size();i++)
      for(int t(0);t<props[i].size();t++)
	{
	  int tt ;
	  double re, im ;
	  file>>tt;
	  for(int l(0);l<props[i][t].size();l++){
	    file>>re; //props[i][t][l].re ;
	    if(fu==0)
	      props[i][t][l] = DComplex(re,0.0) ;
	    else{
	      file>>im ; //props[i][t][l].im ;
	      props[i][t][l]=DComplex(re,im) ;
	    }
	  }
	  if(t!=tt){
	    cerr<<"Wrong data in file: "<<f.c_str() ;
	    throw ;
	  }
	}
    
    file.close();

    return Nx ;
  }


  void WriteProplist(Array<Array< Array<double> > >& props, 
		     int Nx, const string& f){
    ofstream file(f.c_str());
    // Write header line                                                  
    int nprops(props.size());
    int Nt(props[0].size());
    file<<nprops <<" "<<Nt<<" 0 "<<Nx<<" "<<props[0][0].size()<<endl;
    for(int i(0);i<props.size();i++)
      for(int t(0);t<props[i].size();t++){
        file<<t<<" ";
	for(int l(0);l<props[i][t].size();l++)
	  file<<props[i][t][l]<<" ";
	file<<endl ;
      }
    file.close();
  }

  void WriteProplist(Array<Array< Array<DComplex> > >& props, 
		     int Nx, const string& f){
    ofstream file(f.c_str());
    // Write header line                                                  
    int nprops(props.size());
    int Nt(props[0].size());
    file<<nprops <<" "<<Nt<<" 1 "<<Nx<<" "<<props[0][0].size()<<endl;
    for(int i(0);i<props.size();i++)
      for(int t(0);t<props[i].size();t++){
        file<<t<<" ";
	for(int l(0);l<props[i][t].size();l++)
	  file<<props[i][t][l].real()<<" "<<props[i][t][l].imag()<<" " ;
	file<<endl ;
      }
    file.close();
  }

  int ReadProplist(Array< Array< Array2d<double> > >& pr, const std::string& file){
    Array< Array< Array<double> > > tt ;
    int r  = ReadProplist(tt, file) ;
    //pr = tt ;
    pr.resize(tt.size());
    for (int i(0);i<pr.size();i++){
      pr[i].resize(tt[i].size());
      for (int k(0);k<pr[i].size();k++){
	int ss = int(sqrt((double)tt[i][k].size())) ;
	if( ss*ss != tt[i][k].size()){
	  std::cerr<<"Array2d.size(): "<<ss<<std::endl;
	  std::cerr<<"Conversion to square matrix failed. ";
	  std::cerr<<"Array.size(): "<<tt[i][k].size()<<std::endl;
	}
	pr[i][k].resize(ss, ss);   // always resize
	for (int c(0);c<ss;c++)
	  for (int d(0);d<ss;d++)
	    pr[i][k][c][d] = tt[i][k][d+ss*c];
      }
    }
      
    return r ;
  }

  void WriteProplist(Array< Array< Array2d<double> > >& pr, int Nx, const std::string& file){
    Array< Array< Array<double> > > tt ;
    //tt=pr ;  // the underlying equal operators are somehow messed up...
    std::cerr<<"void WriteProplist(Array< Array< Array2d<double> > >& pr, int Nx, const std::string& file): does not work\n" ;
    exit(1);
    WriteProplist(tt, Nx, file) ;
  }

  int ReadProplist(Array< Array< Array2d<DComplex> > >& pr, const std::string& file){
    Array< Array< Array<DComplex> > > tt ;
    int r  = ReadProplist(tt, file) ;
    //pr = tt ;
    pr.resize(tt.size());
    for (int i(0);i<pr.size();i++){
      pr[i].resize(tt[i].size());
      for (int k(0);k<pr[i].size();k++){
	int ss = int(sqrt((double)tt[i][k].size())) ;
	if( ss*ss != tt[i][k].size()){
	  std::cerr<<"Array2d.size(): "<<ss<<std::endl;
	  std::cerr<<"Conversion to square matrix failed. ";
	  std::cerr<<"Array.size(): "<<tt[i][k].size()<<std::endl;
	} 
	pr[i][k].resize(ss, ss);   // always resize
	for (int c(0);c<ss;c++)
	  for (int d(0);d<ss;d++)
	    pr[i][k][c][d] = tt[i][k][d+ss*c];
      }
    }
      
    return r ;
  }



  void WriteProplist(Array< Array< Array2d<DComplex> > >& pr, int Nx, const std::string& file){
    Array< Array< Array<DComplex> > > tt ;
    //tt=pr ; // the underlying equal operators are somehow messed up...
    std::cerr<<"void WriteProplist(Array< Array< Array2d<DComplex> > >& pr, int Nx, const std::string& file): does not work\n" ;
    exit(1);
    WriteProplist(tt, Nx, file) ;
  }

}
