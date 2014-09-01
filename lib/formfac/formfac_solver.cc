// $Id: formfac_solver.cc,v 2.0 2008/12/05 04:43:36 edwards Exp $
//
// Form-factors


#include "ensem/ensem.h"
#include "formfac/formfac_solver.h"
#include "formfac/formfac_data_covar.h"
#include "formfac/formfac_ensem_linfit.h"
#include "recipes/nr.h"
#include "formfac/formfac_chisqq.h"

#include <vector>

//using namespace std;


namespace FF 
{
  // Namespace composition
  using namespace ENSEM;

  //! Linear least squares for real systems
  LLSqResults<EnsemReal> llsq(const LLSqComponent<EnsemReal>& sys, double cov_tol, double solve_tol)
  {
    LLSqResults<EnsemReal> ret;

    // Sanity checks

    // Get all the row indices.
    std::vector<LLSqRow_t> rowList = sys.rowList();
    int N = rowList.size();

//    std::cout << "N=" << N << std::endl;
    if (N <= 0)
    {
      std::cerr << __func__ << ": N is nonsense" << std::endl;
      exit(1);
    }

    // Find the names the of unknowns
    const std::vector<std::string> unknowns = sys.unknowns();
    int M = sys.nUnknowns();
//    std::cout << "M=" << M << std::endl;
    if (M <= 0)
    {
      std::cerr << __func__ << ": M is nonsense" << std::endl;
      exit(1);
    }
    if (M > N)
    {
      std::cerr << __func__ << ": warning M > N, so underconstrained" << std::endl;
    }

    //
    // Loop over indices and construct matrix and rhs
    // Keep a copy of ensemble objects.
    //
    // don't include a row if the entire thing is zero
    //
    std::vector< std::vector<EnsemReal> >  mat;
    std::vector<EnsemReal>  rhs;

//    std::cout << "construct mat" << std::endl;
    for(int n=0; n < N; ++n)
    {
      const LLSqRow_t& row = rowList[n];

      if (sys.isRowZero(row))
	continue;

      EnsemReal rhs_temp = sys.rhs(row);
      rhs.push_back(rhs_temp);
      
      std::vector<EnsemReal> mat_row;
      for(int m=0; m < M; ++m)
      {
	EnsemReal mat_tmp = sys.mat(row,unknowns[m]);
	mat_row.push_back(mat_tmp);
      }
      
      mat.push_back(mat_row);
    } // for n

    int N_nonzero = rhs.size();
    std::cout << "number of non-zero correlator timeslices imported = " << N_nonzero << " from a possible " << N << std::endl;


    // CONSTRUCT A DATA COVARIANCE - a big, dumb thing
    // Use the unscaled data
    std::cout << "construct a covariance of dimension " << N_nonzero << std::endl;
    Array2d<Real> cov = data_covar(rhs);
    std::cout << "constructed a covariance of dimension " << N_nonzero << std::endl;

    std::cout << "invert the covariance" << std::endl;


    // Array2d<Real> inv_cov = data_covar_inv(cov, cov_tol);
    data_covar_inv_t inv = data_covar_inv(cov, cov_tol);
    std::cout << "inverted the covariance" << std::endl;

    // Build data
    Array<EnsemCorrLinGenFitData_t> data(rhs.size());
    for(int n=0; n < rhs.size(); ++n)
    {
      data[n].f = rhs[n];

      data[n].X.resize(M);

      for(int m=0; m < M; ++m)
	data[n].X[m] = mat[n][m];
    }

    //
    // Big linear system call
    //
    EnsemCorrLinGenFitResults_t res = corrLinGenFit(data, inv.inv_cov, solve_tol);
    if (res.coeffs.size() != M)
    {
      std::cerr << __func__ << ": coeffs.size() != M" << std::endl;
      exit(1);
    }

    ret.count_singular_val = res.count_singular_val;

    //change degrees of freedom according to how many singular values were reset in the covraince and in the linear system
    //sing values in the linear system are *nbins
    int sv_linsys = int( double(ret.count_singular_val) / double(sys.nbins())  + 0.5 );
    std::cout << "avg num of sing values = " << sv_linsys << std::endl;

    int ndof = res.ndof - inv.n_reset_sv + sv_linsys;
    ret.ndof  = ndof;

    ret.x     = res.coeffs;
    ret.x_cov = res.coeff_cov;

    ret.chisq = res.chisq;

    //recompute the Q with the modified ndof
    ret.Q = chisqq_zero(ndof, toDouble(res.chisq));
    //    ret.Q     = res.Q;

    return ret;
  }







  //NOT FIXED UP FOR COVARIANCE

#if 0
  //! Linear least squares for complex systems
  LLSqResults llsq(EnsemVectorComplex& x, EnsemReal& chisq, const LLSqComponent<EnsemComplex>& sys)
  {
    LLSqResults ret;

    double tol = 1.0e-6;

    // Sanity checks

    // Find size of RHS
    int N = 1;
    for(int i=0; i < sys.size().size(); ++i)
      N *= sys.size()[i];

//    std::cout << "N=" << N << std::endl;
    if (N == 0)
    {
      ret.success = false;
      return ret;
    }

    // Find number of unknowns
    int M = sys.nUnknowns();
//    std::cout << "M=" << M << std::endl;
    if (M == 0)
    {
      ret.success = false;
      return ret;
    }

    // Output
    EnsemVectorComplex x_tmp;
    EnsemReal chisq_tmp;

    x_tmp.resize(sys.nbins());
    x_tmp.resizeObs(M);
    chisq_tmp.resize(sys.nbins());

    x_tmp = zero;
    chisq_tmp = zero;

    //
    // Loop over indices and construct matrix and rhs
    // Keep a copy of ensemble objects.
    //
    // This requires some thought if this is the best way
    //
    Array< Array<EnsemComplex> >  mat(N);
    Array<EnsemComplex>  rhs(N);

//    std::cout << "construct mat" << std::endl;
    for(int n=0; n < N; ++n)
    {
      Array<int> ind = crtesn(n, sys.size());

      rhs[n] = rescaleEnsemDown(sys.rhs(ind));

      mat[n].resize(M);
      for(int m=0; m < M; ++m)
	mat[n][m] = rescaleEnsemDown(sys.mat(ind,m));

    } // for n


#if 0
    for(int n=0; n < N; ++n)
    {
      char lin[100];
      sprintf(lin, "rhs_%d", n);
      write(std::string(lin), rhs[n]);

      for(int m=0; m < M; ++m)
      {
	sprintf(lin, "mat_%d_%d", n,m);
	write(std::string(lin), mat[n][m]);
      }
    }
#endif


    //
    // Big loop over ensemble bins
    // 
    // QUESTION: do I need to rescale and expose fluctuations here??
    //
    ret.count_singular_val = 0;

    for(int ibin=0; ibin < sys.nbins(); ++ibin)
    {
//      std::cout << "ibin=" << ibin << std::endl;

      // The objects need by SVD
      NR::Mat_IO_DP  a(2*N,2*M);   // a copy for chisq since "u" is destroyed
      NR::Mat_IO_DP  u(2*N,2*M);
      NR::Vec_IO_DP  b(2*N);
      NR::Vec_IO_DP  xx(2*M);

      // A major hassle - need to fake a real matrix from complex
      // Loop over indices and construct matrix and rhs
      for(int n=0; n < N; ++n)
      {
	b[2*n]   = toDouble(real(peekEnsem(rhs[n], ibin)));
	b[2*n+1] = toDouble(imag(peekEnsem(rhs[n], ibin)));
	  
	for(int m=0; m < M; ++m)
	{
	  a[2*n][2*m]     = toDouble(real(peekEnsem(mat[n][m], ibin)));
	  a[2*n+1][2*m]   = toDouble(imag(peekEnsem(mat[n][m], ibin)));
	  a[2*n][2*m+1]   = -a[2*n+1][2*m];
	  a[2*n+1][2*m+1] = a[2*n][2*m];
	}
      } // for n

      // keep a copy
      for(int n=0; n < a.nrows(); ++n)
	for(int m=0; m < a.ncols(); ++m)
	  u[n][m] = a[n][m];


      // Call SVD
      NR::Vec_IO_DP  w(2*M);
      NR::Mat_IO_DP  v(2*M,2*M);
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
//	  std::cerr << "Found a singular value w=" << w[m] << " at ibin=" << ibin << std::endl;
	  ++ret.count_singular_val;
	  w[m] = 0.0;
	}

      // Back-substitute
      svbksb(u,w,v,b,xx);

      // Evaluate chisq
      double chi = 0.0;
      for(int n=0; n < a.nrows(); ++n)
      {
	double sum = 0.0;
	for(int m=0; m < a.ncols(); ++m)
	  sum += a[n][m] * xx[m];

	double tmp = b[n] - sum;
	chi += tmp*tmp;
      } // for n

      pokeEnsem(chisq_tmp, Real(chi), ibin);

      // Pull out the solution vectors
      VectorComplex tmp;
      tmp.resizeObs(M);

      for(int m=0; m < M; ++m)
	pokeObs(tmp, cmplx(Real(xx[2*m]),Real(xx[2*m+1])), m);

      pokeEnsem(x_tmp, tmp, ibin);
    } // for ibin

    // Rescale fluctuations back
    x = rescaleEnsemUp(x_tmp);
    chisq = rescaleEnsemUp(chisq_tmp);

    return ret;
  }
#endif

} // namespace FF
