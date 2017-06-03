// Spectrum program reading directly from databases (db) 
//

#include <iostream>
#include <sstream>
#include <cstdio>
#include "AllConfStoreDB.h"

#include <covfit/statistics.h>
#include <covfit/Function.h>
#include <covfit/fitters/fitfunctions.h>
#include <covfit/fitter.h>
#include <covfit/FitParams_io.h>

using namespace std;
using namespace CovFit ;
using namespace FILEDB ;

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

int main(int argc, char **argv)
{

  if (argc != 2)
  {
    cerr << "Usage: " << argv[0] << " <params>" ; 
    cerr << endl;
    exit(1);
  }

  PropList corr ;

  DBFitParams_t param ;
  
  
  XMLReader xml_in(argv[1]);
  try{
    read(xml_in,"/param_DB", param) ;
  }
  catch(const string& error) {
    cerr<<"OOOPS!! : "<<error<<"\n" ;
    exit(22);
  }
  xml_in.close();
  
  cout<<"Number of states: "<<param.sta.size()<<endl ;
  

  //initialize fitlist files
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
  for(int s(0);s<param.sta.size();s++){
    Array<ofstream> ftl(param.sta[s].fit.mass_param.size());
    for(int i(0);i<param.sta[s].fit.mass_param.size();i++){
      stringstream fname ;
      fname<<"fitlist-"<<i;
      ftl[i].open(fname.str().c_str(),ofstream::app) ;
      ftl[i]<<"\n"<<param.sta[s].name<<" FITS"<<endl ;
    }
    cout<<"Reading state: "<<param.sta[s].name
	<<" from file: "<<param.sta[s].filename<<endl ;

    string dbase = param.sta[s].filename ;

    //AllConfStoreDB< SerialDBKey<KeyHadron2PtCorr_t>,  SerialDBData<PropList> > db;
    ConfDataStoreDB< SerialDBKey<KeyHadron2PtCorr_t>,  SerialDBData<PropList> > db;
    try
      {
	// Open DB

	if (db.open(dbase, O_RDONLY, 0400) != 0)
	  {
	    std::cerr << __func__ << ": error opening dbase= " << dbase << std::endl;
	    exit(1);
	  }
      }
    catch(const std::string& e) 
      {
	std::cerr << __func__ << ": Caught Exception: " << e << std::endl;
	exit(1);
      }
    catch(std::exception& e) 
      {
	std::cerr << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
	exit(1);
      }

    //    ReadProplist(corr,param.sta[s].filename);
    SerialDBKey<KeyHadron2PtCorr_t> k ;
    SerialDBData<PropList> v ;
    k.key()=param.sta[s].key ;
    std::cout<<" Reading propagator key="<<k.key()<<std::endl ;
    db.get(k,v);
    std::cout<<" Got the data v.data().size()"<<v.data().size()<<std::endl ;
    corr=v.data();
    db.close();
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
    Array<double> time(Nt) ;
    for(int t(0);t<Nt;t++) time[t] = t ;
    CovarMat cm(Nt,param.cov.block);
    cm.SetDataList(time,corr);
    cm.CalcCovarMat();
    
    cm.WriteCovarMat("covar."+param.sta[s].name);
    
    
    /** DEBUG
	for(int i(0);i<param.sta[s].range.size();i++)
	cout<<param.sta[s].range[i]<<endl ;
	XMLFileWriter xml_out("foo.xml");
	write(xml_out,"param_db",param);
	xml_out.close();
	cout<<param.sta[s].const_params.size()<<endl;
    **/
    
    // CONTINUE with the fitting code
    coFitter efit(*expo,e_args);
    Array<double> guesspar, fitpar, fiterr ;
    guesspar = param.sta[s].fit.fit_params ;
    //string tt("fitlist."+param.sta[s].name) ;
    //ofstream ftl(tt.c_str());
    for(int i(0);i<param.sta[s].range.size();i++){
      int mind(param.sta[s].fit.min_dist) ;
      int maxd(param.sta[s].fit.min_dist+param.sta[s].range[i]) ;
      for(;maxd<=param.sta[s].stopper;mind++,maxd++){
	efit.setRange(mind,maxd) ;
	efit.setParams(guesspar);
	
	cout<<"\n\n****************\n" ;  
	efit.SetCovarMat(cm) ;
	cout<<efit.Ndof()<<" PARAMETERS: DISTANCE "<<mind<<" to "<<maxd<<endl ;
	efit.Fit() ;
	efit.getParams(fitpar,fiterr);
	// Need to reorder fitpar and fiterr
	orderMasses(fitpar, fiterr, param.sta[s].fit.mass_param);
	// ... so that the smaller mass appears first

	for(int f(0);f<param.sta[s].fit.mass_param.size();f++){
	  ftl[f]<<mind<<" "<<maxd<<"\t";
	  ftl[f].precision(5);
	  ftl[f]<<scientific<<fitpar[param.sta[s].fit.mass_param[f]]<<"\t" ;
	  ftl[f].precision(2);
	  ftl[f]<<fiterr[param.sta[s].fit.mass_param[f]]<<"\t";
	  
	  ftl[f]<<scientific<<efit.ChiSq()<<"\t";
	  ftl[f]<<efit.Ndof()<<"\t" ;
	  ftl[f]<<fixed<<efit.Confidence() ;
	  if(!efit.Convergence()) ftl[f]<<" *" ;
	  if(efit.FailedFit()) ftl[f]<<" !" ;
	  ftl[f]<<endl ;
	  ftl[f].precision(5);
	}
	ostringstream ss ;
	ss<<"fit_graph-"<<param.sta[s].name<<"-"<<mind<<"-"<<maxd<<".ax" ;
	char tt[100] ;
	strcpy(tt,ss.str().c_str());
	efit.graph(tt);
      }
    }

    for(int i(0);i<param.sta[s].fit.mass_param.size();i++)
      ftl[i].close();
    
    delete expo ;
  }


}
