#include <cmath>
#include "nr.h"
using namespace std;

//void usrfun(NR::Vec_I_DP &x, NR::Vec_O_DP &fvec, NR::Mat_O_DP &fjac);

namespace NR {


  /*  void mnewt(const int ntrial, Vec_IO_DP &x, const DP tolx, const DP tolf)
  {
    int k,i;
    DP errx,errf,d;
    
    int n=x.size();
    Vec_INT indx(n);
    Vec_DP p(n),fvec(n);
    Mat_DP fjac(n,n);
    for (k=0;k<ntrial;k++) {
      usrfun(x,fvec,fjac);
      errf=0.0;
      for (i=0;i<n;i++) errf += fabs(fvec[i]);
      if (errf <= tolf) return;
      for (i=0;i<n;i++) p[i] = -fvec[i];
      ludcmp(fjac,indx,d);
      lubksb(fjac,indx,p);
      errx=0.0;
      for (i=0;i<n;i++) {
	errx += fabs(p[i]);
	x[i] += p[i];
      }
      if (errx <= tolx) return;
    }
    return;
  };
  */

  void mnewt( const ScalFuncMultiVars &f , const int ntrial,  Vec_IO_DP &x, const DP tolx, const DP tolf){
    
    int k,i;
    DP errx,errf,d;
    
    int n=x.size();
    Vec_INT indx(n);
    Vec_DP p(n),fvec(n);
    Mat_DP fjac(n,n);
    for (k=0;k<ntrial;k++) {
      //  usrfun(x,fvec,fjac);
 
      fvec = f.values(x);
      fjac = f.jacobian(x);

      errf=0.0;
      for (i=0;i<n;i++) errf += fabs(fvec[i]);
      if (errf <= tolf) return;
      for (i=0;i<n;i++) p[i] = -fvec[i];
      ludcmp(fjac,indx,d);
      lubksb(fjac,indx,p);
      errx=0.0;
      for (i=0;i<n;i++) {
	errx += fabs(p[i]);
	x[i] += p[i];
      }
      if (errx <= tolx) return;
    }
    return;


  }; 

}
