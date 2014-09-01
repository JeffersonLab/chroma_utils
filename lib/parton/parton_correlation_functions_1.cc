// $Id: parton_correlation_functions_1.cc,v 2.0 2008/12/05 04:43:37 edwards Exp $
/*! \file
 *  \brief Correlation functions
 *
 * These routines extract the data from a 3-pt Manage object.  All other
 * variants of correlation functions reference the data through only these calls
 * (excluding Raw_Q_C calls).
 */

#include "parton/parton_correlation_functions.h"
#include "formfac/formfac_manage_3pt.h"
#include <assert.h>

#define _DST_ 4

namespace FF 
{

/*#####################################################################################*/
/*#####################################################################################*/

  //! Get a raw 3-pt
  EnsemVectorComplexF Raw_Q_C(Manage3PtFuncReduced& threept,
			      const PiPf& pi_pf,
			      const ArrayInt& links,
			      const int gamma)
  {
#if DEBUG_CorrelationFunctions_1 == 1
    {
      cout << "BEGIN: " << __FILE__ << ": " << __func__ << endl;
      cout << "p_i = " << pi_pf.p_i << "  p_f = " << pi_pf.p_f << std::endl;
      cout << "links = " << c << std::endl;
    }
#endif

    ThreePtArgReduced arg;
    arg.pi_pf = pi_pf;
    arg.gamma = gamma;
    arg.links = links;

    // Snarf the data
    EnsemVectorComplexF thr = threept[arg];

    // some basic sanity checks
    assert( gamma < 16 );
    if ( isNaN(thr) )
    {
      std::cerr << __func__ << ": found a nan" << std::endl;
      exit(1);
    }

#if DEBUG_CorrelationFunctions_1 == 1
    {
      cout << "END: " << __FILE__ << ": " << __func__ << std::endl;
    }
#endif

    return thr;
  }

/*#####################################################################################*/
/*#####################################################################################*/

  //! \f$\bar{q}\Gamma q\f$
  EnsemVectorComplexF QBar_Gamma_Q(Manage3PtFuncReduced& threept,
				   const PiPf& pi_pf,
				   int gamma)
  {
    ArrayInt mu;  // empty
    return Raw_Q_C(threept, pi_pf, mu, gamma);
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorComplexF QBar_Gamma_UDag_Q(Manage3PtFuncReduced& threept,
					const PiPf& pi_pf,
					int gamma,
					int mu1)
  {
    ArrayInt mu(1);
    mu[ 0 ] = mu1;

    return Raw_Q_C(threept, pi_pf, mu, gamma);
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorComplexF QBar_Gamma_UDagUDag_Q(Manage3PtFuncReduced& threept,
					    const PiPf& pi_pf,
					    int gamma,
					    int mu2,
					    int mu1)
  {
    ArrayInt mu(2);

    /* If mu1 and mu2 cancel, then collapse the resulting pair of U and UDag into just 1. */
    if( abs( mu1 - mu2 ) == _DST_ )
    {
      return QBar_Gamma_Q(threept, pi_pf, gamma);
    }

    mu[ 0 ] = mu1;  mu[ 1 ] = mu2;

    return Raw_Q_C(threept, pi_pf, mu, gamma);
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorComplexF QBar_Gamma_UDagUDagUDag_Q(Manage3PtFuncReduced& threept,
						const PiPf& pi_pf,
						int gamma,
						int mu3,
						int mu2,
						int mu1)
  {
    ArrayInt mu(3);

    /* If mu1 and mu2 cancel, then collapse the resulting pair of U and UDag into just 1. */
    if( abs( mu1 - mu2 ) == _DST_ )
    {
      return QBar_Gamma_UDag_Q(threept, pi_pf, gamma, mu3 );
    }

    /* If mu2 and mu3 cancel, then collapse the resulting pair of U and UDag into just 1. */
    if( abs( mu2 - mu3 ) == _DST_ )
    {
      return QBar_Gamma_UDag_Q(threept, pi_pf, gamma, mu1);
    }

    mu[ 0 ] = mu1;  mu[ 1 ] = mu2;  mu[ 2 ] = mu3;

    return Raw_Q_C(threept, pi_pf, mu, gamma);
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorComplexF QBar_Gamma_UDagUDagUDagUDag_Q(Manage3PtFuncReduced& threept,
						    const PiPf& pi_pf,
						    int gamma,
						    int mu4,
						    int mu3,
						    int mu2,
						    int mu1)
  {
    /* produces warnings */
    /*const direction mu[ 4 ] = { mu1, mu2, mu3, mu4 };*/
    /*const unsigned long int index = LinkPatternIndex( 4, mu );*/

    ArrayInt mu(4);
    mu[ 0 ] = mu1;  mu[ 1 ] = mu2;  mu[ 2 ] = mu3;  mu[ 3 ] = mu4;

    /* check for pairs of mu's which cancel */
    std::cerr << __func__ << ": Four Link Objects: Not Yet Implemented" << std::endl;
    exit(1);

    return Raw_Q_C(threept, pi_pf, mu, gamma);
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorComplexF QBar_Gamma_UDagUDagUDagUDagUDag_Q(Manage3PtFuncReduced& threept,
							const PiPf& pi_pf,
							int gamma,
							int mu5,
							int mu4,
							int mu3,
							int mu2,
							int mu1)
  {
    /* produces warnings */
    /*const direction mu[ 5 ] = { mu1, mu2, mu3, mu4, mu5 };*/
    /*const unsigned long int index = LinkPatternIndex( 5, mu );*/

    ArrayInt mu(5);
    mu[ 0 ] = mu1;  mu[ 1 ] = mu2;  mu[ 2 ] = mu3;  mu[ 3 ] = mu4;  mu[ 4 ] = mu5;

    /* check for pairs of mu's which cancel */
    std::cerr << __func__ << ": Five Link Objects: Not Yet Implemented" << std::endl;
    exit(1);

    return Raw_Q_C(threept, pi_pf, mu, gamma);
  }

/*#####################################################################################*/
/*#####################################################################################*/

} // namespace FF
