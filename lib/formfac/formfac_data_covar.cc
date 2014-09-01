// $Id: formfac_data_covar.cc,v 2.0 2008/12/05 04:43:35 edwards Exp $

#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "formfac/formfac_data_covar.h"
#include "recipes/nr.h"

//#include <ieeefp.h>
#include <assert.h>

using namespace ENSEM;

namespace FF
{

  //-----------------------------------------------------------------------------------
  //! Covariance on real data
  Array2d<Real> data_covar(const std::vector<EnsemReal>& data)
  {
    const int N = data.size();
    Array2d<Real> covar(N,N);

    // The means
    std::vector<Real> meas(N);
    for(int i=0; i < N; ++i)
      meas[i] = mean(data[i]);

    // The covariance itself
    bool uncorrP = false;
    const int nbins = data[0].size();

    for(int i=0; i < N; ++i)
    {
      for(int j=i; j < N; ++j)
      {
	Real sum(0.0);
	for(int ibin=0; ibin < nbins; ++ibin)
	{
	  Real d_i = peekEnsem(data[i], ibin);
	  Real d_j = peekEnsem(data[j], ibin);

	  sum += d_i*d_j;
	}

	sum        /= Real(nbins);
	covar(i,j)  = (sum - meas[i]*meas[j]) / Real(nbins - 1);

	if (uncorrP && (i != j))
	  covar(i,j) = 0.0;

	covar(j,i) = covar(i,j);
      }
    }

    return covar;
  } // void func 


  //-----------------------------------------------------------------------------------
  //! Invariance of a subset of the covariance matrix
  /*! */
  //Array2d<Real> data_covar_inv(const Array2d<Real>& covar, double tol, int Tmin, int Tmax)
  data_covar_inv_t data_covar_inv(const Array2d<Real>& covar, double tol, int Tmin, int Tmax)
  {
    int N = covar.size1();
    if (covar.size2() != N)
    {
      std::cerr << __func__ << ": covar not square" << std::endl;
      exit(1);
    }
    
    Array2d<Real> inv_covar(N,N);

    assert(Tmin >= 0); assert(Tmin < N);
    assert(Tmax >= 0); assert(Tmax < N);
    assert(Tmin <= Tmax);

    // The objects need by SVD
    NR::Mat_IO_DP  u(N,N);
    NR::Vec_IO_DP  sig(N);

    // Initialize
    for(int i=0; i < N; ++i)
      for(int j=0; j < N; ++j)
      {
	u[i][j] = 0.0;
      }

    for(int i=0; i < N; ++i)
    {
      u[i][i] = 1.0;
      sig[i]  = 1.0;
    }

    // Rescale by diagonal elements to stabilize
    for(int n=Tmin; n <= Tmax; ++n)
      sig[n] = toDouble(sqrt(covar(n,n)));


    // The inverse of the whole matrix. However, only a subset
    // is non-trivial
    // Loop over indices and construct the subset of the matrix
    for(int n=Tmin; n <= Tmax; ++n)
    {
      for(int m=Tmin; m <= Tmax; ++m)
      {
	u[n][m] = toDouble(covar(n,m)) / (sig[n]*sig[m]);
      }
    } // for n


    // Call SVD
    NR::Vec_IO_DP  w(N);
    NR::Vec_IO_DP  b(N);
    NR::Vec_IO_DP  xx(N);
    NR::Mat_IO_DP  v(N,N);

    svdcmp(u, w, v);

    // Find and reset singular values
    NR::DP  wmax = 0.0;
    for(int m=Tmin; m <= Tmax; ++m)
    {
      if (w[m] > wmax) 
	wmax = w[m];
    }

    int count_singular_val = 0;
    NR::DP  wmin = wmax * tol;
    for(int m=Tmin; m <= Tmax; ++m)
    {
//    std::cout << __func__ << ": w[" << m << "]= " << w[m] << std::endl;
      if (w[m] < wmin)
      {
//	std::cerr << __func__ << ": Found a singular value w= " << w[m] << std::endl;
	++count_singular_val;
	w[m] = 0.0;
      }
    }

    if (count_singular_val > 0)
      std::cout << __func__ << ": Number of singular values= " << count_singular_val << std::endl;

    // Back-substitute for each column.
    // Yes, I know. It can be problematic to compute the entire
    // matrix inverse rather than back-substitute whenever it is needed.
    // However, structurally, it's much cleaner to pass around a matrix
    for(int m=0; m < N; ++m)
    {
      // Initialize the column unit-vector
      for(int n=0; n < N; ++n)
      {
	b[n] = 0.0;
      }

      b[m] = 1.0;

      // Solve for the column vector
      svbksb(u,w,v,b,xx);
    
      // Evaluate chisq
      for(int n=0; n < N; ++n)
      {
	inv_covar(n,m) = Real(xx[n]);
      }
    } // for(m)

    // Scale back out the rescaling
    for(int n=Tmin; n <= Tmax; ++n)
    {
      for(int m=Tmin; m <= Tmax; ++m)
      {
	inv_covar(n,m) /= Real(sig[n]*sig[m]);
      }
    } // for n

    //build the output struct
    data_covar_inv_t out;
    out.inv_cov = inv_covar;
    out.n_reset_sv = count_singular_val;
    out.Tmin = Tmin;
    out.Tmax = Tmax;
    //and set it free to frolic in land of linear systems 
    return out;
  } 

 //-----------------------------------------------------------------------------------
  //! Invert the whole of the covariance matrix
  /*! */
  data_covar_inv_t data_covar_inv(const Array2d<Real>& covar, double tol)
  {
    return data_covar_inv(covar, tol, 0, covar.size1() - 1);
  }








  //-----------------------------------------------------------------------------------
  //! Write correlation matrix
  void writeCorrMat(XMLWriter& xml_out, const std::string& path, const Array2d<Real>& covar)
  {
    push(xml_out, path);
    for(int i=0; i < covar.size1(); ++i)
    {
      push(xml_out, "elem");
      write(xml_out, "i", i);

      for(int j=0; j < covar.size2(); ++j)
      {
	push(xml_out, "elem");
	write(xml_out, "j", j);
	write(xml_out, "covar", covar(i,j));
	pop(xml_out); // elem
      }
      
      pop(xml_out); // elem
    }
    pop(xml_out); // Covariance
  }

 


  //-----------------------------------------------------------------------------------
  //! Write correlation matrix
  void writeCorrMatTest(XMLWriter& xml_out, const std::string& path, 
			const Array2d<Real>& covar,
			const Array2d<Real>& inv_covar,
			int Tmin, int Tmax)
  {
    const int N= covar.size1();
    Array2d<Real> ident(N,N);

    // Init
    for(int i=0; i < N; ++i)
      for(int j=0; j < N; ++j)
	ident(i,j) = zero;

    // Multiply covar*inv_covar
    for(int i=Tmin; i <= Tmax; ++i)
      for(int j=Tmin; j <= Tmax; ++j)
	for(int k=Tmin; k <= Tmax; ++k)
	  ident(i,j) += covar(i,k)*inv_covar(k,j);

    // Print
    push(xml_out, path);
    for(int i=0; i < covar.size1(); ++i)
    {
      push(xml_out, "elem");
      write(xml_out, "i", i);

      for(int j=0; j < covar.size2(); ++j)
      {
	push(xml_out, "elem");
	write(xml_out, "j", j);
	write(xml_out, "ident", ident(i,j));
	pop(xml_out); // elem
      }
      
      pop(xml_out); // elem
    }
    pop(xml_out); // Covariance
  }


} // end namespace
