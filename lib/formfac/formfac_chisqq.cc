// $Id: formfac_chisqq.cc,v 2.0 2008/12/05 04:43:35 edwards Exp $
/*! \file
 * \brief Chisq
 */

#include "recipes/nr.h"

using namespace std;

namespace FF
{

  // ********************** FUNCTION CHISQQ *****************************
  //
  // Returns (logically)  gammq(ndof/2,chisq/2)
  //
  double chisqq(int n, double x)
  {
    if(x < 0 || n <= 0)
    {
      cerr << __func__ << ": Error in CHISQQ! Need n.gt.0 and x.ge.0 !";
      cerr << " n= " << n << "  x= " << x << endl;
      exit(1);
    }

    return ( (double(n)*log(x)-x) < 2.0*log(1.0e-10) ) ? 
      NR::gammq(double(n)*0.50, x*0.5) :
      1.0-NR::gammp(double(n)*0.50, x*0.5);
  }


  //
  // Returns (logically)  gammq(ndof/2,chisq/2)
  //
  double chisqq_zero(int ndof, double x)
  {
    return ( (ndof > 0) && ( x >=0 ) ) ? chisqq(ndof, x) : 0.0;
  }
} // end namespace

