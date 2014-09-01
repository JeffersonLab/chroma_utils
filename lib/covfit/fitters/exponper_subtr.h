//$Id: exponper_subtr.h,v 2.0 2008/12/05 04:43:33 edwards Exp $
//$Log: exponper_subtr.h,v $
//Revision 2.0  2008/12/05 04:43:33  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:19  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.1  2007/05/10 02:16:06  kostas
//added the fit function I need for periodic two meson spectroscopy
//
/*! Definition of the sum of exponentials with periodic BC function */


#ifndef EXPONENTIALS_PERIODIC_SUBTR
#define EXPONENTIALS_PERIODIC_SUBTR
//#define FUNCTION_NAME Exponentials
#include "covfit/Function.h"

namespace CovFit {

class ExponentialsPeriodicSubtr: public Function{
private:
  double center ; // the t=0 point
  double period ; // periodic extent of the lattice
public:
  
  ExponentialsPeriodicSubtr():Function(){}
  ExponentialsPeriodicSubtr(int n):center(0.0),period(0.0),Function(2*n){f_init();}
  ExponentialsPeriodicSubtr(double c):center(c),period(0.0),Function(){f_init();}
  ExponentialsPeriodicSubtr(double p, int n):center(0.0),period(p),Function(2*n){f_init();}
  ExponentialsPeriodicSubtr(double c, double p):center(c),period(p),Function(){f_init();}
  ExponentialsPeriodicSubtr(double c, double p,int n):center(c),period(p),Function(2*n){f_init();}
  
  virtual ~ExponentialsPeriodicSubtr(){}

  void setCenter(double c){
    center = c ;
  }

  void f_init(){
    if(npar%2!=0)
      std::cerr<<"ExponentialsPeriodicSubtr::f_init : npar has to be even\n";
  }
  void f_init(FILE *fp)
  {
    int n ;
    fprintf(stderr,"sum of exponentials fitting with periodic b.c.\n");
    fprintf(stderr,"The value in the midle of the lattice is subtracted.\n");
    fprintf(stderr,"This is needed for two pion states.\n");
    fprintf(stderr,"Enter the center:\n");
    fscanf(fp,"%lg",&center);
    fprintf(stderr,"Enter the periodic extent of the lattice:\n");
    fscanf(fp,"%lg",&period);
    fprintf(stderr,"Enter the number of exponentials:\n");
    fscanf(fp,"%d",&n);
    npar = 2*n ;
    fprintf(stderr,"Using center: %g\n",center);
    fprintf(stderr,"Number of parameters: %i\n",npar);
  }
  
  double f(double x, double *b)
  {
    double z;
    double xx(fabs(x-center) - period/2.0);
    double yy(period-xx);
    int i;
    z=0;
    for(i=0;i<npar;i+=2) z += b[i]*(exp(-b[i+1]*xx)+exp(-b[i+1]*yy) - 2.0);
    return(z);
  }
  
  double df(double x, double *b, int i)
  {
    double xx(fabs(x-center)-period/2.0);
    double yy(period-xx);
    if(i%2)
      return -b[i-1]*(xx*exp(-b[i]*xx)+yy*exp(-b[i]*yy));
    else 
      return exp(-b[i+1]*xx)+exp(-b[i+1]*yy) -2.0 ;
  }
  
  double ddf(double x, double *b, int i, int j)
  {
    double xx(fabs(x-center)-period/2.0);
    double yy(period-xx);
    if(i%2)
      if(j==i)
	return b[i-1]*(xx*xx*exp(-b[i]*xx)+yy*yy*exp(-b[i]*yy)) ;
      else if(j==i-1)
	return -xx*exp(-b[i]*xx)-yy*exp(-b[i]*yy) ;
      else
	return 0.0 ;
    else
      if(j==i+1)
	return  -xx*exp(-b[i+1]*xx)-yy*exp(-b[i+1]*yy) ;
      else
	return 0.0 ;
  }
  
};

} // namespace CovFit

#endif
