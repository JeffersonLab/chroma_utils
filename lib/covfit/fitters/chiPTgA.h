//$Id: chiPTgA.h,v 2.0 2008/12/05 04:43:32 edwards Exp $
//$Log: chiPTgA.h,v $
//Revision 2.0  2008/12/05 04:43:32  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:19  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.2  2006/12/06 04:17:58  kostas
//added function "factory" ... Robert sorry for the abuse of terminology...
//
//Revision 1.1  2006/07/12 00:20:26  kostas
//added some things about gA fits.
//
//Revision 1.1  2005/11/23 03:18:51  kostas
//moved the fitter header files into the main include tree
//
//Revision 1.3  2005/10/05 14:11:26  kostas
//added mpi^4 term in chiPTfpi and also added it to the makefiles
//
//Revision 1.2  2005/10/03 16:04:02  kostas
//fixed bugs
//
//Revision 1.1  2005/10/03 15:52:52  kostas
//fpi chipt fit function
//
/*! ChiPT formula for gA
*/


#ifndef CHIPT_gA
#define CHIPT_gA
#include "covfit/Function.h"

namespace CovFit {

class ChiPTgA: public Function{
private:
  gAJfunc J0 ;
  gAJfunc J  ;
  gAJfunc K  ;

  double _PI_ ;
  double inv_PI_SQ ;
  

public:
  
  ChiPTgA():Function(3, double D):J0(0.0,1.0),J(D,1.0),K(D,1.0){
    npar=5 ; 

    _PI_= 3.14159265358979323848 ;
    inv_PI_SQ = 1.0/(_PI_*_PI_);
  }


  

  virtual ~ChiPTgA(){}


  void f_init(){
    
  }
  void f_init(FILE *fp)
  {
    fprintf(stderr,"One loop chiPT for gA\n");
    fprintf(stderr,"gA = gA_0 +... Beane and Savage \n");
    fprintf(stderr,"b[0] = gA_0   b[1] = gND  b[2] = gDD  b[3] = fpi \n b[4]=C the counterterm\n");
  }

#define gA_0 b[0] 
#define gND  b[1]
#define gDD  b[2]    
#define fpi  b[3]
#define CNN  b[4] 
#define gDD2  (b[2]*b[2])    
#define fpi2  (b[3]*b[3])



  // x is mpi^2 b[0] = f b[1]= quadratic term 
  double f(double x, double *b)
    {
      double x2 = x*x ;
      
      double gA = b[0] - 0.5*inv_PI_SQ/fpi2*gA_0*gA_0*gA_0*J0.f(x2,NULL) ;
      gA -=  0.5*inv_PI_SQ/fpi2*gND2*(gA_0 + 25.0/81.0*gDD)*J.f(x2,NULL) ;
      gA -=  0.25*inv_PI_SQ/fpi2*gA_0*x2*(log(x2) - 1.0) ;
      gA +=  4.0/9.0*inv_PI_SQ/fpi2*gND2*gA_0*K.f(x2,NULL) + CNN*x2;
      return gA ;
    }
  
  //derivatives not done yet.
  double df(double x, double *b, int i)
    {
      if(i==0)
	return 1.0 + inv_EIGHT_PI_SQ*x/(b[0]*b[0])*log(x) ;
      else if( i==1 )
	return x ;
      else
	return x*x ;
    }
  
  double ddf(double x, double *b, int i, int j)
  {
    if((i==0) && (j == 0))
      return - 2.0 * inv_EIGHT_PI_SQ*x/(b[0]*b[0]*b[0])*log(x) ;
    else
      return 0.0 ;
  }
  
};

} // namespace CovFit

#endif
