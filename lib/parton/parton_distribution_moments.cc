// $Id: parton_distribution_moments.cc,v 2.0 2008/12/05 04:43:38 edwards Exp $
/*! \file
 *  \brief Parton distribution moments
 */

#include "parton/parton_distribution_moments.h"
#include "parton/parton_derivatives.h"
#include <assert.h>

namespace FF
{

/*#####################################################################################*/
/* gamma matrix conventions                                                            */
/*#####################################################################################*/

  // anonymous namespace
  namespace
  {
    const RealD twopi(6.28318530717958647688);

    const int dX = 0;
    const int dY = 1;
    const int dZ = 2;
    const int dT = 3;

    const int i_I        =  0;
    const int i_GX       =  1;
    const int i_GY       =  2;
    const int i_GXGY     =  3;
    const int i_GZ       =  4;
    const int i_GXGZ     =  5;
    const int i_GYGZ     =  6;
    const int i_GXGYGZ   =  7;
    const int i_GT       =  8;
    const int i_GXGT     =  9;
    const int i_GYGT     = 10;
    const int i_GXGYGT   = 11;
    const int i_GZGT     = 12;
    const int i_GXGZGT   = 13;
    const int i_GYGZGT   = 14;
    const int i_GXGYGZGT = 15;
    
    const int i_mG5GXGYGZGT =  0;
    const int i_G5GYGZGT    =  1;
    const int i_mG5GXGZGT   =  2;
    const int i_G5GZGT      =  3;
    const int i_G5GXGYGT    =  4;
    const int i_mG5GYGT     =  5;
    const int i_G5GXGT      =  6;
    const int i_mG5GT       =  7;
    const int i_mG5GXGYGZ   =  8;
    const int i_G5GYGZ      =  9;
    const int i_mG5GXGZ     = 10;
    const int i_G5GZ        = 11;
    const int i_G5GXGY      = 12;
    const int i_mG5GY       = 13;
    const int i_G5GX        = 14;
    const int i_mG5I        = 15;
  } // namespace anonymous

//**************************************************************************************
// This factor stuff seems funky and out of place. For now, set the mass
// and energy to one. That's how Dru initialized them in his extract code.

  // anonymous namespace
  namespace
  {
    const RealF aE(1.0);
    const RealF aM(1.0);
  }


/*#####################################################################################*/
/*#####################################################################################*/

  Array<RealF> vec_pf(Manage3PtFuncReduced& threept, const PiPf& pi_pf)
  {
    Array<RealF> P_f(pi_pf.p_f.size());

    for(int mu=0,n=0; mu < threept.lattSize().size(); ++mu)
    {
      if (mu != threept.decayDir())
      {
	P_f[n] = twopi * RealF(pi_pf.p_f[n]) / RealF(threept.lattSize()[mu]);
	++n;
      }
    }

    return P_f;
  }

/*#####################################################################################*/
/*#####################################################################################*/
/*                                                                                     */
/* The prefactor is P_t / S_z which for P_z = 0 and spin in the z direction gives      */
/* E / M = 1.  For P_x = P_y = 0 we have P_t = gamma * M = E and S_z = gamma * M = E,  */
/* which gives E / E.  For now we assume that P_z = 0.                                 */
/*                                                                                     */
/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorRealF q_d1(Manage3PtFuncReduced& threept, const PiPf& pi_pf)
  {
    assert( pi_pf.p_f[2] == 0 );
      
    EnsemVectorRealF D;
    
    D  = imag(QBar_Gamma_D_Q(threept, pi_pf, i_G5GZ, dT)) + imag(QBar_Gamma_D_Q(threept, pi_pf, i_mG5GT, dZ));
    D *= RealF(-2.0) / aM;

    return D;
  }
  
/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorRealF q_d2(Manage3PtFuncReduced& threept, const PiPf& pi_pf)
  {
    std::cerr << __func__ << ": not implemented" << std::endl;
    exit(1);

    EnsemVectorRealF D;
    return D;
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorRealF q_lpol_0(Manage3PtFuncReduced& threept, const PiPf& pi_pf)
  {
    assert( pi_pf.p_f[2] == 0 );

    EnsemVectorRealF D;
    
    D  = imag(QBar_Gamma_Q(threept, pi_pf, i_G5GZ));
    D *= aE / aM;

    return D;
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorRealF q_lpol_1_a(Manage3PtFuncReduced& threept, const PiPf& pi_pf)
  {
    assert( pi_pf.p_f[2] == 0 );
    assert( pi_pf.p_f[0] != 0 );

    Array<RealF> P_f = vec_pf(threept, pi_pf);
    EnsemVectorRealF D;

    D  = real(QBar_Gamma_D_Q(threept, pi_pf, i_G5GX, dZ)) + real(QBar_Gamma_D_Q(threept, pi_pf, i_G5GZ, dX));
    D *= - aE / aM / P_f[0];

    return D;
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorRealF q_lpol_1_b(Manage3PtFuncReduced& threept, const PiPf& pi_pf)
  {
    assert( pi_pf.p_f[2] == 0 );

    EnsemVectorRealF D;

    D  = imag(QBar_Gamma_D_Q(threept, pi_pf, i_G5GZ, dT)) - imag(QBar_Gamma_D_Q(threept, pi_pf, i_mG5GT, dZ));
    D *= RealF(-1) / aM;

    return D;
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorRealF q_lpol_2(Manage3PtFuncReduced& threept, const PiPf& pi_pf)
  {
    assert( pi_pf.p_f[2] == 0 );

    Array<RealF> P_f = vec_pf(threept, pi_pf);
    EnsemVectorRealF D;

    D  = real(QBar_Gamma_DD_Q(threept, pi_pf, i_G5GX, dZ, dT)) +
         real(QBar_Gamma_DD_Q(threept, pi_pf, i_G5GX, dT, dZ)) +
         real(QBar_Gamma_DD_Q(threept, pi_pf, i_G5GZ, dX, dT)) +
         real(QBar_Gamma_DD_Q(threept, pi_pf, i_G5GZ, dT, dX)) -
         real(QBar_Gamma_DD_Q(threept, pi_pf, i_mG5GT, dX, dZ)) -
         real(QBar_Gamma_DD_Q(threept, pi_pf, i_mG5GT, dZ, dX));

    D *= RealF(1.0) / RealF(2.0) / aM / P_f[0];

    return D;
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorRealF q_tpol_0(Manage3PtFuncReduced& threept, const PiPf& pi_pf)
  {
    assert( pi_pf.p_f[2] == 0 );

    return imag(QBar_Gamma_Q(threept, pi_pf, i_G5GZGT));
  }

/*#####################################################################################*/
/*#####################################################################################*/

/* Check the factor of 1/2 again. I think it is right, but Dolgov didn't have it. */

  EnsemVectorRealF q_tpol_1(Manage3PtFuncReduced& threept, const PiPf& pi_pf)
  {
    assert( pi_pf.p_f[2] == 0 );
    assert( pi_pf.p_f[0] != 0 );

    Array<RealF> P_f = vec_pf(threept, pi_pf);
    EnsemVectorRealF D;

    D  = real(QBar_Gamma_D_Q(threept, pi_pf, i_G5GZGT, dX)) + real(QBar_Gamma_D_Q(threept, pi_pf, i_mG5GXGZ, dT));
    D *= - RealF(1.0) / RealF(2.0) / P_f[0];

    return D;
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorRealF q_upol_0(Manage3PtFuncReduced& threept, const PiPf& pi_pf)
  {
    return real(QBar_Gamma_Q(threept, pi_pf, i_GT));
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorRealF q_upol_1_a(Manage3PtFuncReduced& threept, const PiPf& pi_pf)
  {
    assert( pi_pf.p_f[0] != 0 );

    Array<RealF> P_f = vec_pf(threept, pi_pf);
    EnsemVectorRealF D;

    D  = imag(QBar_Gamma_D_Q(threept, pi_pf, i_GX, dT)) + imag(QBar_Gamma_D_Q(threept, pi_pf, i_GT, dX));
    D *= RealF(1) / RealF(2) / P_f[0];

    return D;
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorRealF q_upol_1_b(Manage3PtFuncReduced& threept, const PiPf& pi_pf)
  {
    Array<RealF> P_f = vec_pf(threept, pi_pf);
    EnsemVectorRealF D;

    D  = real(QBar_Gamma_D_Q(threept, pi_pf, i_GT, dT)) - 
      (  real(QBar_Gamma_D_Q(threept, pi_pf, i_GX, dX)) +
	 real(QBar_Gamma_D_Q(threept, pi_pf, i_GY, dY)) +
	 real(QBar_Gamma_D_Q(threept, pi_pf, i_GZ, dZ)) ) / RealF(3.0);

    D *= -aE / ( aE*aE - RealF(1)/RealF(3) *( P_f[0]*P_f[0] + P_f[1]*P_f[1] + P_f[2]*P_f[2] ) );

    return D;
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorRealF q_upol_2(Manage3PtFuncReduced& threept, const PiPf& pi_pf)
  {
    Array<RealF> P_f = vec_pf(threept, pi_pf);
    EnsemVectorRealF D;

    D  =
      real(QBar_Gamma_DD_Q(threept, pi_pf, i_GX, dX, dT)) +
      real(QBar_Gamma_DD_Q(threept, pi_pf, i_GX, dT, dX)) +
      real(QBar_Gamma_DD_Q(threept, pi_pf, i_GT, dX, dX)) - 
      RealF(0.5) * ( real(QBar_Gamma_DD_Q(threept, pi_pf, i_GY, dY, dT)) +
		     real(QBar_Gamma_DD_Q(threept, pi_pf, i_GY, dT, dY)) +
		     real(QBar_Gamma_DD_Q(threept, pi_pf, i_GT, dY, dY)) +
		     real(QBar_Gamma_DD_Q(threept, pi_pf, i_GZ, dZ, dT)) +
		     real(QBar_Gamma_DD_Q(threept, pi_pf, i_GZ, dT, dZ)) +
		     real(QBar_Gamma_DD_Q(threept, pi_pf, i_GT, dZ, dZ)) );

    D *= -RealF(1) / RealF(3) / ( P_f[0]*P_f[0] - RealF(0.5) * ( P_f[1]*P_f[1] + P_f[2]*P_f[2] ) );

    return D;
  }

/*#####################################################################################*/
/*#####################################################################################*/

  EnsemVectorRealF q_upol_3(Manage3PtFuncReduced& threept, const PiPf& pi_pf)
  {
    std::cerr << __func__ << ": not implemented" << std::endl;
    exit(1);
    EnsemVectorRealF D;
    return D;
  }

} // namespace FF
