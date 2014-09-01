// -*- C++ -*-
// $Id: adat_lapack.h,v 2.0 2008/12/05 04:43:32 edwards Exp $


#ifndef ADAT_LAPACK_H
#define ADAT_LAPACK_H
//#ifdef MAC_OSX
//MAC OS X
//#include <vecLib/vBLAS.h>
//#include <veclib/clapack.h>
//#endif

#include <covfit/proplist.h> 
extern "C" {
  /** symmetic generalized eigenvalue problem **/
int dsygv_( int    *ITYPE, 
	    char   *JOBZ, 
	    char   *UPLO, 
	    int    *N, 
	    Double *A, 
	    int    *LDA, 
	    Double *B, 
	    int    *LDB, 
	    Double *W, 
	    Double *WORK,
	    int    *LWORK, 
	    int    *INFO ) ;
  
}
#endif
