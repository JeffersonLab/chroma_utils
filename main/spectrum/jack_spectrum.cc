// $Id: jack_spectrum.cc,v 2.0 2008/12/05 04:44:02 edwards Exp $
// $Log: jack_spectrum.cc,v $
// Revision 2.0  2008/12/05 04:44:02  edwards
// Changed to version 2.0.
//
// Revision 1.2  2007/04/16 21:24:57  kostas
// added jack_spectrum
//
// Revision 1.1  2007/04/15 04:22:53  kostas
// jackknife lists from spectrum fits
//
//

#include <iostream>
#include <sstream>
#include <cstdio>

#include <covfit/statistics.h>
#include <covfit/Function.h>
#include <covfit/fitters/fitfunctions.h>
#include <covfit/fitter.h>
#include <covfit/FitParams_io.h>

using namespace std;
using namespace CovFit ;

//assuming that these are exponentials so we reorder the
//amplitudes together with the masses
void orderMasses( Array<Double>& ofp,  
		  Array<Double>& oefp,
		  const Array<int>& im){
  for(int j(0);j<im.size()-1;j++){
    for(int i(0);i<im.size()-j-1;i++){
      if(ofp[im[i]]>ofp[im[i+1]]){//Reorder
        Double A = ofp[im[i]-1] ;
        Double M = ofp[im[i]  ] ;
	Double eA = oefp[im[i]-1] ;
        Double eM = oefp[im[i]  ] ;

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

int main(int argc, char **argv)
{

  if (argc != 2)
  {
    cerr << "Usage: " << argv[0] << " <params>" ; 
    cerr << endl;
    exit(1);
  }

  PropList corr ;

  FitParams_t param ;
  
  
  XMLReader xml_in(argv[1]);
  try{
    read(xml_in,"/param_db", param) ;
  }
  catch(const string& error) {
    cerr<<"OOOPS!! : "<<error<<"\n" ;
    exit(22);
  }
  xml_in.close();
  
  cout<<"Number of states: "<<param.sta.size()<<endl ;

  //initialize fitlist files
  /**
  for(int s(0);s<param.sta.size();s++){
    Array<ofstream> ftl(param.sta[s].fit.mass_param.size());
    for(int i(0);i<param.sta[s].fit.mass_param.size();i++){
      stringstream fname ;
      fname<<"fitlist-"<<i;
      ftl[i].open(fname.str().c_str()) ;
      ftl[i]<<endl;
      ftl[i].close();
    }
  }
  **/

  for(int s(0);s<param.sta.size();s++){

    cout<<"Reading state: "<<param.sta[s].name
	<<" from file: "<<param.sta[s].filename<<endl ;

    int Nx = ReadProplist(corr,param.sta[s].filename);
    
    int Ncnfs(corr.size()) ;
    int Nt(corr[0].size());
    
    
    Function* expo;
    try{
      expo = CreateFunction(param.sta[s].fit.fitfunc, 
			    param.sta[s].fit.fit_params.size());
    }
    catch(string e){
      cout<<"OOOPS! "<<e<<endl ;
    }
    try{
      expo = CreateFunction(param.sta[s].fit.fitfunc, param.sta[s].fit.period, 
			    param.sta[s].fit.fit_params.size());
    }
    catch(string e){
      cout<<"OOOPS! "<<e<<endl ;
    }


    //Exponentials expo(param.sta[s].fit_params.size()/2);
    FitterArgs e_args ;
    
    
    //e_args.xlow = (4) ;
    //e_args.xhigh = (17);
    e_args.maxit = param.fit.Niter ;
    e_args.eps   = param.fit.Toler ;
    e_args.npar = expo->Npar();
    e_args.fixpar = param.sta[s].fit.fixedpar ;
    
    //calculate the covariance matrix
    Array<Double> time(Nt) ;
    for(int t(0);t<Nt;t++) time[t] = t ;
    CovarMat cm(Nt,param.cov.block);
    cm.SetDataList(time,corr);
    cm.CalcCovarMat();
    
    cm.WriteCovarMat("covar."+param.sta[s].name);
    
    
    
    // CONTINUE with the fitting code
    coFitter efit(*expo,e_args);
    Array<Double> guesspar, fitpar, fiterr ;
    guesspar = param.sta[s].fit.fit_params ;
    //string tt("fitlist."+param.sta[s].name) ;
    //ofstream ftl(tt.c_str());
    //for(int i(0);i<param.sta[s].range.size();i++){
      int mind(param.sta[s].fit.min_dist) ;
      int maxd(param.sta[s].fit.max_dist) ;
      
      efit.setRange(mind,maxd) ;

      Array< Array<Double> > jfitpar(Ncnfs) ;
      for(int j(0);j<Ncnfs;j++)
	{

	  PropList jcorr ;
	  splice(jcorr,corr,j,j+param.cov.block-1) ;

	  cm.SetDataList(time,jcorr);
	  cm.CalcCovarMat();

	  efit.setParams(guesspar);
	
	  cout<<"\n\n****************\n" ;  
	  efit.SetCovarMat(cm) ;
	  cout<<efit.Ndof();
	  cout<<" PARAMETERS: DISTANCE "<<mind<<" to "<<maxd<<endl ;
	  efit.Fit() ;
	  efit.getParams(jfitpar[j],fiterr);
	  // Need to reorder fitpar and fiterr
	  orderMasses(jfitpar[j], fiterr, param.sta[s].fit.mass_param);
	  // ... so that the smaller mass appears first
	  guesspar = jfitpar[j];
	}//jackknife loop
      
      //RESCALE THE FLUCTUATIONS
      for(int i(0);i<param.sta[s].fit.fit_params.size();i++){
	PropList  par(Ncnfs);
	stringstream fname ;
	fname<<"fitpar-"<<i<<"_"<<param.sta[s].name<<".jk";	
	for(int j(0);j<Ncnfs;j++){
	  par[j].resize(1);
	  par[j][0] = jfitpar[j][i] ;	
	}
	Array<Double> mpar = mean(par) ;
	//Now do the rescaling
	for(int j(0);j<Ncnfs;j++){
	  Double d(jfitpar[j][i] - mpar[0]) ;
	  par[j][0] = mpar[0] + Double(Ncnfs-1)*d ;
	}
	string tt(fname.str().c_str());
	WriteProplist(par,Nx,tt) ;
      }
      
      delete expo ;
  }


}
