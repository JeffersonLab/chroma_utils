//$Id: chiPT_fKofPi.h,v 2.0 2008/12/05 04:43:32 edwards Exp $
//$Log: chiPT_fKofPi.h,v $
//Revision 2.0  2008/12/05 04:43:32  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:19  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.1  2008/07/03 15:48:02  kostas
//*** empty log message ***
//
//
/*! ChiPT formula for fK/fPi
*/


#ifndef CHIPT_mKomPi
#define CHIPT_mKomPi
#include "covfit/Function.h"
#include "Lookup.h"

namespace CovFit {

class ChiPTfKofPi: public Function{
private:
  double inv_EIGHT_PI_SQ ;
  LookUp mKoFpi ;

public:
  
  ChiPTfKofPi(const Lookup& t):Function(3),mKoFpi(t){
    npar=3 ; 
    inv_SIXTEEN_PI_SQ = 1.0/157.91367041742973790310;

  }
  
  virtual ~ChiPTfKofPi(){}


  void f_init(){
    
  }
  void f_init(FILE *fp)
  {
    fprintf(stderr,"One loop chiPT for mK/mPi\n");
    fprintf(stderr,"fK/fPi  =") ;
    fprintf(stderr,"1 + logs + e2 x + e4 x^2\n");
    fprintf(stderr,"b[0] = e2   b[1] = e4  \n");
  }
  
  // x is mpi^2 b[0] = f b[1]= quadratic term 
  double f(double x, double *b)
    {
      double Mk = mKoFpi.f(x,NULL) ;
      double Mk2 = Mk*Mk ;

      double meta2 = 4.0/3.0*mK2 - 1.0/3.0*x ;
      double mu_pi = inv_SIXTEEN_PI_SQ*x*log(x) ;
      double mu_eta =  inv_SIXTEEN_PI_SQ*meta2*log(meta2) ;
      double mu_K =  inv_SIXTEEN_PI_SQ*Mk2*log(Mk2) ;

      return ( 1 + 5.0/4.0*mu_pi - 0.5*mu_K - 3.0/4.0*mu_eta + 
	       b[0]*(Mk2 - x) +  b[1]*(Mk2 - x)*(Mk2 - x)) ;

    }
  
  double df(double x, double *b, int i)
    {
      double Mk = mKoFpi.f(x,NULL) ;
      double Mk2 = Mk*Mk ;
      if(i==0)
	return (Mk2 - x) ;
      else
	return (Mk2 - x)*(Mk2 - x)  ;
    }
  
  double ddf(double x, double *b, int i, int j){
    return 0.0 ;
  }
  
};

} // namespace CovFit

#endif
