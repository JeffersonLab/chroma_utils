// $Id: formfac_ensem_linfit.cc,v 2.0 2008/12/05 04:43:35 edwards Exp $
/*! \file
 * \brief Chi-sq fit to a sum of linear params with generic function weights
 */

#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "formfac/formfac_ensem_linfit.h"
#include "formfac/formfac_linear_covar_genfit.h"
#include "formfac/formfac_chisqq.h"
#include <assert.h>

using namespace ENSEM;

namespace FF
{

  //-----------------------------------------------------------------------------------
  //  Fit to the single value function  f(\vec{x}) = \sum_i a_i * X_i(\vec{x}) with covariance
  EnsemCorrLinGenFitResults_t corrLinGenFit(const Array<EnsemCorrLinGenFitData_t>& data,
					    const Array2d<Real>& inv_covar,
					    double tol,
					    int Tmin, int Tmax)
  {
    EnsemCorrLinGenFitResults_t res;

    assert(Tmin >= 0); assert(Tmin < data.size());
    assert(Tmax >= 0); assert(Tmax < data.size());
    assert(Tmin <= Tmax);

    if (inv_covar.size1() != inv_covar.size2())
    {
      std::cerr << __func__ << ": inv_covar is not square" << std::endl;
      exit(1);
    }
    if (data[Tmin].X.size() <= 0)
    {
      std::cerr << __func__ << ": no fit functions" << std::endl;
      exit(1);
    }

    // Relevant sizes
    const int N     = data[Tmin].X.size();
                            //std::cout << "number of form-factors is " << N << std::endl;
    const int nbins = data[Tmin].f.size();
    //  std::cout << "data has "  << nbins << " bins" << std::endl;

    for(int n=0; n < N; ++n)
    {
      //  std::cout << "basis has " << data[Tmin].X[n].size() << " bins" << std::endl;
      if (data[Tmin].X[n].size() != nbins)
      {
	std::cerr << __func__ << ": fit function ensembles not consistent size" << std::endl;
	exit(1);
      }
    }

    // Make a copy of the data but in rescaled format
    Array<EnsemCorrLinGenFitData_t> data_resc(data.size());
    for(int i=Tmin; i <= Tmax; ++i)
    {
      data_resc[i].f = rescaleEnsemDown(data[i].f);

      data_resc[i].X.resize(N);
      for(int n=0; n < N; ++n)
	data_resc[i].X[n] = rescaleEnsemDown(data[i].X[n]);
    }

    //
    // Jackknife loops
    //
    res.coeffs.resize(N);
    for(int n=0; n < N; ++n)
      res.coeffs[n].resize(nbins);
    
    res.coeff_cov.resize(N,N);
    for(int n=0; n < N; ++n)
      for(int m=0; m < N; ++m)
	res.coeff_cov[n][m].resize(nbins);

    res.chisq = zero;
    res.count_singular_val = 0;


    for(int ibin=0; ibin < nbins; ++ibin)
    {
      Array<CorrLinGenFitData_t> data_bin(data_resc.size());
      for(int i=Tmin; i <= Tmax; ++i)
      {
	data_bin[i].f = peekEnsem(data_resc[i].f, ibin);

	data_bin[i].X.resize(N);
	for(int n=0; n < N; ++n)
	  data_bin[i].X[n] = peekEnsem(data_resc[i].X[n], ibin);
      }
      
      CorrLinGenFitResults_t res_bin = corrLinGenFit(data_bin, inv_covar, tol, Tmin, Tmax);

      res.ndof   = res_bin.ndof;
      res.chisq += res_bin.chisq;
      res.count_singular_val += res_bin.count_singular_val;

      for(int n=0; n < N; ++n)
	{
	  pokeEnsem(res.coeffs[n], res_bin.coeffs[n].a, ibin);
	  for(int m=0; m < N; ++m){pokeEnsem(res.coeff_cov[n][m], res_bin.coeff_cov[n][m], ibin);}
	}

    } // for(ibin)

    for(int n=0; n < N; ++n)
      {
	res.coeffs[n] = rescaleEnsemUp(res.coeffs[n]);
	for(int m=0; m < N; ++m){ res.coeff_cov[n][m] = rescaleEnsemUp(res.coeff_cov[n][m]);}
      }
    // And the confidence level
    res.chisq /= Real(nbins);
    res.Q = chisqq_zero(res.ndof, toDouble(res.chisq));

    return res;
  }


  //  Fit to the single value function  f(\vec{x}) = \sum_i a_i * X_i(\vec{x}) with covariance
  EnsemCorrLinGenFitResults_t corrLinGenFit(const Array<EnsemCorrLinGenFitData_t>& data,
					    const Array2d<Real>& inv_covar,
					    double tol)
  {
    return corrLinGenFit(data, inv_covar, tol, 0, data.size()-1);
  }

} // end namespace
