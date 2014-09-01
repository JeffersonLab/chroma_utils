// $Id: parton_derivatives.cc,v 2.0 2008/12/05 04:43:38 edwards Exp $
/*! \file
 *  \brief Derivative constructions
 */

#include "parton/parton_derivatives.h"
#include "parton/parton_correlation_functions.h"

#define _DST_ 4

namespace FF
{

  namespace
  {
    const Double twopi(6.28318530717958647688);
  }


/*#####################################################################################*/
/*#####################################################################################*/

  Array<ComplexF> expiqmu(Manage3PtFuncReduced& threept, const PiPf& pi_pf)
  {
    Array<ComplexF> eiqmu(4);
    Array<int> q = pi_pf.p_f - pi_pf.p_i;

    for(int mu=0,n=0; mu < eiqmu.size(); ++mu)
    {
      if (mu != threept.decayDir())
      {
	RealF theta = twopi * RealF(q[n++]) / RealF(threept.lattSize()[mu]);
	eiqmu[mu] = cmplx(cos(theta), sin(theta));
      }
      else
      {
	eiqmu[mu] = cmplx(RealF(1), RealF(0));
      }
    }

    return eiqmu;
  }


/*#####################################################################################*/
/*#####################################################################################*/

  Array<ComplexF> expmiqmu(Manage3PtFuncReduced& threept, const PiPf& pi_pf)
  {
    Array<ComplexF> emiqmu(4);
    Array<ComplexF> eiqmu = expiqmu(threept, pi_pf);

    for(int mu=0; mu < eiqmu.size(); ++mu)
    {
      emiqmu = conj(eiqmu[mu]);
    }

    return emiqmu;
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorComplexF QBar_Gamma_D_Q(Manage3PtFuncReduced& threept,
				     const PiPf& pi_pf,
				     int gamma,
				     int mu1)
  {
    const int delta_t_mu1 = ( mu1 == 3 ) ? 1 : 0;
    const int minus_mu1 = mu1 + _DST_;

    /* these are the various phase factors required to achieve a momentum projection at */
    /* the center of the derivative */
    Array<ComplexF>  eiqmu = expiqmu(threept, pi_pf);
    Array<ComplexF> emiqmu = expmiqmu(threept, pi_pf);

    EnsemVectorComplexF D;

    /* frwd mu1 */
    D  = emiqmu[mu1] * cshift(QBar_Gamma_UDag_Q(threept, pi_pf, gamma, minus_mu1), + delta_t_mu1);
    D -=  eiqmu[mu1] * cshift(QBar_Gamma_UDag_Q(threept, pi_pf, gamma, mu1),       - delta_t_mu1);

    /* bkwd mu1 */
    D -= QBar_Gamma_UDag_Q(threept, pi_pf, gamma, mu1);
    D += QBar_Gamma_UDag_Q(threept, pi_pf, gamma, minus_mu1);

    /* normalize the 4 terms */
    D *= RealF(1) / RealF(4);

    return D;
  }


/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorComplexF QBar_Gamma_DD_Q(Manage3PtFuncReduced& threept,
				      const PiPf& pi_pf,
				      int gamma,
				      int mu2,
				      int mu1)
  {
    const int delta_t_mu1 = ( mu1 == 3 ) ? 1 : 0;
    const int delta_t_mu2 = ( mu2 == 3 ) ? 1 : 0;
    const int minus_mu1 = mu1 + _DST_;
    const int minus_mu2 = mu2 + _DST_;

    /* these are the various phase factors required to achieve a momentum projection at */
    /* the center of the derivative */
    Array<ComplexF>  eiqmu = expiqmu(threept, pi_pf);
    Array<ComplexF> emiqmu = expmiqmu(threept, pi_pf);

    EnsemVectorComplexF D;

    /* frwd mu2 - frwd mu1 */
    D  = emiqmu[mu2] * emiqmu[mu1] * cshift(QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma, 
								  minus_mu2, minus_mu1),
					    + delta_t_mu2 + delta_t_mu1);
    D -= emiqmu[mu2] *  eiqmu[mu1] * cshift(QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma, 
								  minus_mu2,       mu1),
					    + delta_t_mu2 - delta_t_mu1);
    D -=  eiqmu[mu2] * emiqmu[mu1] * cshift(QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma, 
								        mu2, minus_mu1), 
					    - delta_t_mu2 + delta_t_mu1);
    D +=  eiqmu[mu2] *  eiqmu[mu1] * cshift(QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma,
								        mu2,       mu1), 
					    - delta_t_mu2 - delta_t_mu1);

    /* frwd mu2 - bkwd mu1 */
    D -= emiqmu[mu2] * cshift(QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma, 
						    minus_mu2,       mu1), 
			      + delta_t_mu2);
    D += emiqmu[mu2] * cshift(QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma, 
						    minus_mu2, minus_mu1), 
			      + delta_t_mu2);
    D +=  eiqmu[mu2] * cshift(QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma,
						    mu2,       mu1),
			      - delta_t_mu2);
    D -=  eiqmu[mu2] * cshift(QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma,
						    mu2, minus_mu1),
			      - delta_t_mu2);

    /* bkwd mu2 - frwd mu1 */
    D -= emiqmu[mu1] * cshift(QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma,
						    mu2, minus_mu1),
			      + delta_t_mu1);
    D +=  eiqmu[mu1] * cshift(QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma,
						    mu2,       mu1), 
			      - delta_t_mu1);
    D += emiqmu[mu1] * cshift(QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma,
						    minus_mu2, minus_mu1), 
			      + delta_t_mu1);
    D -=  eiqmu[mu1] * cshift(QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma,
						    minus_mu2,       mu1),
			      - delta_t_mu1);
  
    /* bkwd mu2 - bkwd mu1 */
    D += QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma,       mu2,       mu1);
    D -= QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma, minus_mu2,       mu1);
    D -= QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma,       mu2, minus_mu1);
    D += QBar_Gamma_UDagUDag_Q(threept, pi_pf, gamma, minus_mu2, minus_mu1);

    /* normalize the 16 terms */
    D *= RealF(1) / RealF(16);

    return D;
  }


/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorComplexF QBar_Gamma_DDD_Q(Manage3PtFuncReduced& threept,
				       const PiPf& pi_pf,
				       int gamma,
				       int mu3,
				       int mu2,
				       int mu1)
  {
    EnsemVectorComplexF D;

    std::cerr << __func__ << ": not implemented yet" << std::endl;
    exit(1);

    return D;
  }


} // namespace FF
