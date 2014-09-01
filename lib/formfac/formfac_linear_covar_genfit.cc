// $Id: formfac_linear_covar_genfit.cc,v 2.0 2008/12/05 04:43:35 edwards Exp $
/*! \file
 * \brief Chi-sq fit to a sum of linear params with generic function weights
 */

#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "formfac/formfac_linear_covar_genfit.h"
#include "formfac/formfac_chisqq.h"

#include "recipes/nr.h"
#include <assert.h>

using namespace ENSEM;
using namespace std;


namespace FF
{

  //  Fit to the single value function  f(\vec{x}) = \sum_i a_i * X_i(\vec{x}) with covariance
  CorrLinGenFitResults_t corrLinGenFit(const Array<CorrLinGenFitData_t>& data,
				       const Array2d<Real>& inv_covar,
				       double tol,
				       int Tmin, int Tmax)
  {
    CorrLinGenFitResults_t res;

    if (inv_covar.size1() != inv_covar.size2())
    {
      std::cerr << __func__ << ": inv_covar is not square" << std::endl;
      exit(1);
    }
    if (data.size() != inv_covar.size1())
    {
      std::cerr << __func__ << ": data and inv_covar are not the same size" << std::endl;
      exit(1);
    }

    assert(Tmin >= 0); assert(Tmin < data.size());
    assert(Tmax >= 0); assert(Tmax < data.size());
    assert(Tmin <= Tmax);

    res.Tmin = Tmin;
    res.Tmax = Tmax;

    // The size of the coordinate space which determines the number of params to fit
    const int N = data[Tmin].X.size();
    int M = N;

    NR::Mat_IO_DP  u(N,M);
    NR::Vec_IO_DP  b(N);

    // Initialize
    for(int n=0; n < u.nrows(); ++n)
    {
      b[n] = 0.0;
      for(int m=0; m < u.ncols(); ++m)
      {
	u[n][m] = 0.0;
      }
    }

    // Accumulate would-be covariance matrix
    for(int i=Tmin; i <= Tmax; ++i)
    {
      const CorrLinGenFitData_t& d_i = data[i];

      for(int j=Tmin; j <= Tmax; ++j)
      {
	const CorrLinGenFitData_t& d_j = data[j];

	// Build the matrix and accumulate
	for(int n=0; n < u.nrows(); ++n)
	{
	  b[n] += toDouble(d_i.X[n] * inv_covar(i,j) * d_j.f);
	  
	  for(int m=0; m < u.ncols(); ++m)
	  {
	    u[n][m] += toDouble(d_i.X[n] * inv_covar(i,j) * d_j.X[m]);
	  }
	}
      }
    }

    // The objects need by SVD
    NR::Vec_IO_DP  xx(M);

    res.count_singular_val = 0;

    // Call SVD
    NR::Vec_IO_DP  w(M);
    NR::Mat_IO_DP  v(M,M);
    svdcmp(u, w, v);

    // Find and reset singular values
    NR::DP  wmax = 0.0;
    for(int m=0; m < w.size(); ++m)
      if (w[m] > wmax) 
	wmax = w[m];

    NR::DP  wmin = wmax * tol;
    for(int m=0; m < w.size(); ++m)
      if (w[m] < wmin)
      {
//	std::cerr << __func__ << ": Found a singular value w=" << w[m] << std::endl;
	++res.count_singular_val;
	w[m] = 0.0;
      }

//    if (res.count_singular_val > 0)
//      std::cout << __func__ << ": Found number singular values= " << res.count_singular_val << std::endl;

    // Back-substitute
    svbksb(u,w,v,b,xx);

    // Get the covariances
    NR::Mat_IO_DP  cvm(M,M);
    svdvar(v, w, cvm);

    // Extract results
    res.coeffs.resize(N);
    for(int n=0; n < N; ++n)
    {
      res.coeffs[n].a     = xx[n];
      res.coeffs[n].a_err = sqrt(cvm[n][n]);
    }

    //record the off-diag cov too
    res.coeff_cov.resize(N, N);
    for(int n=0; n<N; ++n)
       for(int m=0; m<N; ++m)
	 {
	   res.coeff_cov[n][m] = cvm[n][m]; 
	 }
	 

    // Construct deviates
    Array<Real> dev(data.size());
    for(int i=Tmin; i <= Tmax; ++i)
    {
      const CorrLinGenFitData_t& d_i = data[i];

      dev[i] = -d_i.f;
      for(int n=0; n < N; ++n)
	dev[i] += res.coeffs[n].a * d_i.X[n];
    }

    // Evaluate chisq
    res.chisq = zero;

    for(int i=Tmin; i <= Tmax; ++i)
    {
      for(int j=Tmin; j <= Tmax; ++j)
      {
	res.chisq += dev[i] * inv_covar(i,j) * dev[j];
      }
    }

    // And the confidence level
    res.ndof = Tmax-Tmin+1 - M;
    res.Q = chisqq_zero(res.ndof, toDouble(res.chisq));
    
    return res;
  }


  //!  Write the results
  void write(XMLWriter& xml_out, const std::string& xml_group, 
	     const CorrLinGenFitResults_t::Coeff_t& param)

  {
    push(xml_out, xml_group);

    write(xml_out, "a", param.a);
    write(xml_out, "a_err", param.a_err);

    pop(xml_out);
  }


  //  Write the results
  void write(XMLWriter& xml_out, const std::string& xml_group, const CorrLinGenFitResults_t& fit)
  {
    push(xml_out, xml_group);
    write(xml_out, "Tmin", fit.Tmin);
    write(xml_out, "Tmax", fit.Tmax);
    write(xml_out, "coeffs", fit.coeffs);
    write(xml_out, "count_singular_val", fit.count_singular_val);
    write(xml_out, "ndof", fit.ndof);
    write(xml_out, "Q", fit.Q);
    write(xml_out, "chisq", fit.chisq);
    pop(xml_out);
  }


  //  Fit to the single value function  f(\vec{x}) = \sum_i a_i * X_i(\vec{x}) with covariance
  CorrLinGenFitResults_t corrLinGenFit(const Array<CorrLinGenFitData_t>& data,
				       const Array2d<Real>& inv_covar,
				       double tol)
  {
    return corrLinGenFit(data, inv_covar, tol, 0, data.size()-1);
  }


  //  Fit to the single value function  f(\vec{x}) = \sum_i a_i * X_i(\vec{x}) with covariance
  CorrLinGenFitResults_t corrLinGenFit(XMLWriter& xml_out, const std::string& xml_group,
				       const Array<CorrLinGenFitData_t>& data,
				       const Array2d<Real>& inv_covar,
				       double tol)
  {
    CorrLinGenFitResults_t  fit = corrLinGenFit(data, inv_covar, tol);
    write(xml_out, xml_group, fit);
    return fit;
  }


  //  Fit to the single value function  f(\vec{x}) = \sum_i a_i * X_i(\vec{x}) with covariance
  CorrLinGenFitResults_t corrLinGenFit(XMLWriter& xml_out, const std::string& xml_group,
				       const Array<CorrLinGenFitData_t>& data,
				       const Array2d<Real>& inv_covar,
				       double tol,
				       int Tmin, int Tmax)
  {
    CorrLinGenFitResults_t  fit = corrLinGenFit(data, inv_covar, tol, Tmin, Tmax);
    write(xml_out, xml_group, fit);
    return fit;
  }

} // end namespace

