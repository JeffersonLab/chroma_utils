// $Id: parton_correlation_functions_2.cc,v 2.0 2008/12/05 04:43:38 edwards Exp $
/*! \file
 *  \brief Correlation functions
 *
 * These routines extract the data from a 3-pt Manage object.  All other
 * variants of correlation functions reference the data through only these calls
 * (excluding Raw_Q_C calls).
 *
 * The correlation functions in this file (CorrelationFunctions_2.c) achieve retrieve
 * values by calls to the correlation functions in the file CorrelationFunctions_1.c.
 */

#include "parton/parton_correlation_functions.h"
#include "formfac/formfac_manage_3pt.h"

/*#####################################################################################*/
/*#####################################################################################*/

#define ND 4
#define DST 4
#define _2DST 8
#define DSTm1 3
#define _2DSTm1 7

/*#####################################################################################*/
/*#####################################################################################*/

namespace FF 
{

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorComplexF QBar_U_Gamma_Q(Manage3PtFuncReduced& threept,
				     const PiPf& pi_pf,
				     int gamma,
				     int mu1)
  {
    int MinusMu1;
    int TPlusMu1 = 0;

    if( mu1 < DST )
    {
      MinusMu1 = mu1 + DST;

      if( mu1 == DSTm1 )
      {
	TPlusMu1 += 1;
      }
    }
    else
    {
      MinusMu1 = mu1 - DST;

      if( mu1 == _2DSTm1 )
      {
	TPlusMu1 -= 1;
      }
    }

    return cshift(QBar_Gamma_UDag_Q(threept, pi_pf, gamma, MinusMu1), TPlusMu1);
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorComplexF QBar_UU_Gamma_Q(Manage3PtFuncReduced& threept,
				      const PiPf& pi_pf,
				      int gamma,
				      int mu1,
				      int mu2)
  {
    int MinusMu1;
    int MinusMu2;
    int TPlusMu1PlusMu2 = 0;

    if( mu1 < DST )
    {
      MinusMu1 = mu1 + DST;

      if( mu1 == DSTm1 )
      {
	TPlusMu1PlusMu2 += 1;
      }
    }
    else
    {
      MinusMu1 = mu1 - DST;

      if( mu1 == _2DSTm1 )
      {
	TPlusMu1PlusMu2 -= 1;
      }
    }

    if( mu2 < DST )
    {
      MinusMu2 = mu2 + DST;

      if( mu2 == DSTm1 )
      {
	TPlusMu1PlusMu2 += 1;
      }
    }
    else
    {
      MinusMu2 = mu2 - DST;

      if( mu2 == _2DSTm1 )
      {
	TPlusMu1PlusMu2 -= 1;
      }
    }

    return cshift(QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma, MinusMu1, MinusMu2), TPlusMu1PlusMu2);
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorComplexF QBar_UUU_Gamma_Q(Manage3PtFuncReduced& threept,
				       const PiPf& pi_pf,
				       int gamma,
				       int mu1,
				       int mu2,
				       int mu3)
  {
    int MinusMu1;
    int MinusMu2;
    int MinusMu3;
    int TPlusMu1PlusMu2PlusMu3 = 0;

    if( mu1 < DST )
    {
      MinusMu1 = mu1 + DST;

      if( mu1 == DSTm1 )
      {
	TPlusMu1PlusMu2PlusMu3 += 1;
      }
    }
    else
    {
      MinusMu1 = mu1 - DST;

      if( mu1 == _2DSTm1 )
      {
	TPlusMu1PlusMu2PlusMu3 -= 1;
      }
    }

    if( mu2 < DST )
    {
      MinusMu2 = mu2 + DST;

      if( mu2 == DSTm1 )
      {
	TPlusMu1PlusMu2PlusMu3 += 1;
      }
    }
    else
    {
      MinusMu2 = mu2 - DST;

      if( mu2 == _2DSTm1 )
      {
	TPlusMu1PlusMu2PlusMu3 -= 1;
      }
    }

    if( mu3 < DST )
    {
      MinusMu3 = mu3 + DST;

      if( mu3 == DSTm1 )
      {
	TPlusMu1PlusMu2PlusMu3 += 1;
      }
    }
    else
    {
      MinusMu3 = mu3 - DST;

      if( mu3 == _2DSTm1 )
      {
	TPlusMu1PlusMu2PlusMu3 -= 1;
      }
    }

    return cshift(QBar_Gamma_UDagUDagUDag_Q(threept, pi_pf, gamma, MinusMu1, MinusMu2, MinusMu3), TPlusMu1PlusMu2PlusMu3);
  }

/*#####################################################################################*/
/*#####################################################################################*/

} // namespace FF  
