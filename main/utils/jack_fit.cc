// Revision 2.0  2008/12/05 04:44:06  edwards
// Changed to version 2.0.
//
// Revision 1.1  2007/10/06 04:04:25  kostas
// added a jack_fit util to make ensemble files of masses.
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

#include <covfit/FitParams_io.h>

#include <covfit/statistics.h>
#include <covfit/Function.h>
#include <covfit/fitters/fitfunctions.h>
#include <covfit/fitter.h>
#include <covfit/jackknife.h>


//#include "jackknife_io.h"
//#include "basic_functions.h"

using namespace std;
using namespace CovFit ;
//assuming that these are exponentials so we reorder the
//amplitudes together with the masses
void orderMasses( Array<double>& ofp,  
                  Array<double>& oefp,
                  const Array<int>& im){
  for(int j(0);j<im.size()-1;j++){
    for(int i(0);i<im.size()-j-1;i++){
      if(ofp[im[i]]>ofp[im[i+1]]){//Reorder
        double A = ofp[im[i]-1] ;
        double M = ofp[im[i]  ] ;
        double eA = oefp[im[i]-1] ;
        double eM = oefp[im[i]  ] ;

        //suffle masses
        ofp[im[i  ]-1] = ofp[im[i+1]-1] ;
        ofp[im[i  ]  ] = ofp[im[i+1]  ] ;
        ofp[im[i+1]-1] = A ;
        ofp[im[i+1]  ] = M ;

        //suffle errors
        oefp[im[i  ]-1] = oefp[im[i+1]-1] ;
        oefp[im[i  ]  ] = oefp[im[i+1]  ] ;
        oefp[im[i+1]-1] = eA ;
        oefp[im[i+1]  ] = eM ;
      }
    }
  }

}


void write(const std::string& t,
           const Array<double>& v,
           const Array<double>& e){
  //Needs size checking                                                        
  for(int i(0);i<v.size();i++)
    cout<<t<<": "<<i<<" "<<v[i]<<" "<<e[i]<<endl ;
}

int main(int argc, char **argv)
{
  if ( ( argc !=2 ) )
  {
    cerr << "Usage: " << argv[0] << " ";
    cerr<< "<xml input>";
    cerr << endl;
    exit(1);
  }

  FitterParam fit ;
  string in_file ;
  string ensem_out ;
  XMLReader xml_in(argv[1]);
  try{
    read(xml_in,"/param_db/fit", fit) ;
    read(xml_in,"/param_db/filename", in_file) ;
    read(xml_in,"/param_db/ensem_out", ensem_out) ;
  }
  catch(const string& error) {
    cerr<<"OOOPS!! : "<<error<<"\n" ;
    exit(22);
  }
  xml_in.close();

  cout<<"ensem_out: "<<ensem_out<<endl ;
  cout<<"fit_function: "<<fit.fitfunc<<endl ;

  PropList corr ;
  cout<<"Reading state from file: "<<in_file<<endl ;
  
  int Nx = ReadProplist(corr,in_file);
  
  int Ncnfs(corr.size()) ;
  int Nt(corr[0].size());
  
    
  Function* expo;
  try{
    expo = CreateFunction(fit.fitfunc, fit.fit_params.size());
    }
  catch(string e){
    cout<<"OOOPS! "<<e<<endl ;
  }
  try{
    expo = CreateFunction(fit.fitfunc, fit.period, fit.fit_params.size());
  }
  catch(string e){
    cout<<"OOOPS! "<<e<<endl ;
  }

  PropList mass(Ncnfs) ;
  FitterArgs e_args ;

  e_args.maxit = 1000 ;
  e_args.eps   = 1.0e-5 ;
  e_args.npar = expo->Npar();
  e_args.fixpar = fit.fixedpar ;
  
  //calculate the covariance matrix
  Array<double> time(Nt) ;
  for(int t(0);t<Nt;t++) time[t] = t ;
  CovarMat cm(Nt,1);

  coFitter efit(*expo,e_args);
  Array<double> guesspar, fitpar, fiterr ;

  guesspar = fit.fit_params ;

  //Jackknife loop 
  for(int j(0);j<Ncnfs;j++){
    PropList jcorr ;
    splice(jcorr,corr,j,j) ;
    cm.SetDataList(time,jcorr);
    cm.CalcCovarMat();
    
    
    Array<double> jfitpar,fiterr ;
    mass[j].resize(Nt);
    for(int mind(fit.min_dist);mind<fit.max_dist;mind++){
      //guesspar = fit.fit_params ;
      efit.setRange(mind,fit.max_dist) ;
      efit.setParams(guesspar);
      efit.SetCovarMat(cm) ;

      cout<<"\n\n****************\n" ;  
      cout<<efit.Ndof();
      cout<<" PARAMETERS: DISTANCE "<<mind<<" to "<<fit.max_dist<<endl ;
      efit.Fit() ;
      efit.getParams(jfitpar,fiterr);
      orderMasses(jfitpar, fiterr, fit.mass_param);
      guesspar = jfitpar ;
      mass[j][mind] = jfitpar[fit.mass_param[0]] ;
    }
    for(int t(0);t<fit.min_dist;t++)
      mass[j][t] =  mass[j][fit.min_dist];//just to keep Proplist Readers happy
    for(int t(fit.max_dist);t<Nt;t++)
      mass[j][t] = mass[j][fit.max_dist-1];//just to keep Proplist Readers happy
    
  }


  Array<double> av_mass = mean(mass) ;

  PropList tt_mass(Ncnfs) ;
  for(int j(0);j<Ncnfs;j++)
    {
      Array<double> d = mass[j] - av_mass ;
      tt_mass[j] = av_mass + double(Ncnfs-1)*d ;
    }
  WriteProplist(tt_mass,Nx,ensem_out);

}
