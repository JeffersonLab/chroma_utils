/*! Binding momentum function */


#ifndef BOUNDMOM
#define BOUNDMOM
//#define FUNCTION_NAME Exponentials
#include "covfit/Function.h"

namespace CovFit {

class BoundMom: public Function{

public:
  
  BoundMom():Function(2){}
  
  virtual ~BoundMom(){}

  void f_init(){
  }
  void f_init(FILE *fp)
  {
    int n ;
    fprintf(stderr,"Binding Mom vs L  fitting\n");
    fprintf(stderr,"   b[0] + b[1]/x * (exp(-b[0]*x) + sqrt(2)*exp(-sqrt(2)*b[0]*x)) \n");
  }
  
  double f(double x, double *b)
  {
    return (b[0] + b[1]/x*(exp(-b[0]*x) + sqrt(2.0)*exp(-b[0]*sqrt(2.0)*x)) );
  }
  
  double df(double x, double *b, int i)
  {
    if(i==1)
      return 1.0/x*(exp(-b[0]*x) + sqrt(2.0)*exp(-b[0]*sqrt(2.0)*x)) ;
    else if(i==0)
      return 1.0 - b[1]*(exp(-b[0]*x) + 2.0*exp(-b[0]*sqrt(2.0)*x)) ;
    else 
      return (0) ;
  }
  
  double ddf(double x, double *b, int i, int j)
  {
    int ii(i),jj(j) ;
    if(i>j){ 
      ii=j ;
      jj=i ;
    }
    
    if(ii==0)
      if(jj==0)
	return b[1]*x*(exp(-b[0]*x) + 2.0*sqrt(2.0)*exp(-b[0]*sqrt(2.0)*x)) ;
      else if(jj==1)
	return -(exp(-b[0]*x) + 2.0*exp(-b[0]*sqrt(2.0)*x)) ;
      else
	return 0.0 ;
    else 
      return 0.0 ;
  }
  
};

} // namespace CovFit

#endif
