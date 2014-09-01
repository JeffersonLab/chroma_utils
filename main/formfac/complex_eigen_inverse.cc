// $Id: complex_eigen_inverse.cc,v 2.0 2008/12/05 04:43:47 edwards Exp $

#include "ensem/ensem.h"
#include "recipes/nr.h"

using namespace ENSEM;
using namespace std;

//! NxN Matrix solver
Array<Real> solveNxN(const Array2d<Real>& a, const Array<Real>& b)
{
  int N = a.size1();
  int M = a.size2();
  if (N != M)
  {
    cerr << __func__ << ": only supports square systems" << endl;
    exit(1);
  }

  if (b.size() != N)
  {
    cerr << __func__ << ": only supports square systems" << endl;
    exit(1);
  }

  // The objects need by SVD
  NR::Mat_IO_DP  u(N,M);
  NR::Vec_IO_DP  xx(M);
  NR::Vec_IO_DP  bb(N);

  // Make a copy
  for(int n=0; n < u.nrows(); ++n)
  {
    bb[n] = toDouble(b[n]);

    for(int m=0; m < u.ncols(); ++m)
    {
      u[n][m] = toDouble(a[n][m]);
    }
  } // for n

  int count_singular_val = 0;
  
  // Call SVD
  NR::Vec_IO_DP  w(M);
  NR::Mat_IO_DP  v(M,M);
  svdcmp(u, w, v);

  // Find and reset singular values
  NR::DP  wmax = 0.0;
  for(int m=0; m < w.size(); ++m)
    if (w[m] > wmax) 
      wmax = w[m];

  double tol = 1.0e-6;
  NR::DP  wmin = wmax * tol;
  for(int m=0; m < w.size(); ++m)
    if (w[m] < wmin)
    {
      std::cerr << __func__ << ": Found a singular value w=" << w[m] << std::endl;
      ++count_singular_val;
      w[m] = 0.0;
    }

  // Back-substitute
  svbksb(u,w,v,bb,xx);

  // Extract results
  Array<Real>  soln(xx.size());
  for(int n=0; n < xx.size(); ++n)
    soln[n] = xx[n];

  return soln;
}



//! 3x3 Matrix solver
Array<EnsemReal> ensemSolveNxN(const Array2d<EnsemReal>& a, const Array<EnsemReal>& b)
{
  // Setup the result
  Array<EnsemReal> soln = b;

  // Make resized copies of matrix and rhs
  Array2d<EnsemReal> aa = a;
  Array<EnsemReal> bb = b;

  for(int n=0; n < a.size2(); ++n)
  {
    bb[n] = rescaleEnsemDown(b[n]);

    for(int m=0; m < a.size1(); ++m)
    {
      aa[n][m] = rescaleEnsemDown(a[n][m]);
    }
  } // for n

  int nboot = b[0].size();

  // Big jackknife loop
  for(int ibin=0; ibin < nboot; ++ibin)
  {
    Array2d<Real> aaa(aa.size2(), aa.size1());
    Array<Real> bbb(bb.size());

    for(int n=0; n < aa.size2(); ++n)
    {
      bbb[n] = peekEnsem(bb[n], ibin);

      for(int m=0; m < a.size1(); ++m)
      {
	aaa[n][m] = peekEnsem(aa[n][m], ibin);
      }
    }

    Array<Real> xxx = solveNxN(aaa,bbb);
    for(int n=0; n < xxx.size(); ++n)
      pokeEnsem(soln[n], xxx[n], ibin);
  }
  for(int n=0; n < soln.size(); ++n)
    soln[n] =  rescaleEnsemUp(soln[n]);
  return soln;
}



// MAIN PROGRAM INCL COMMAND LINE I/O

int main(int argc, char *argv[])
{
  if (argc != 3)
  {
    cerr << "Usage:"
	 << argv[0] << ":  <dimension>  <time>" 
	 << endl;
    exit(1);
  }

  int dim;
  {
    istringstream val(argv[1]);
    val >> dim;
  }

  int t;
  {
    istringstream val(argv[2]);
    val >> t;
  }

  // Read IMAG PART of eigenvectors
  Array2d<EnsemReal> B(dim,dim);
  for(int i=0; i < dim; ++i)
  {
    for(int j=0; j < dim; ++j)
    {
      ostringstream file_name;
      file_name << "ev_" << i << "_" << j << "_im.jknf";

      EnsemVectorReal tmp;
      read(file_name.str(), tmp);
      B[i][j] = peekObs(tmp, t);
    }
  }

  int nbins = B[0][0].size();

  // Invert IMAG PART of eigenvector matrix
  Array2d<EnsemReal> Binv(dim,dim);
  for(int i=0; i < dim; ++i)
  {
    Array<EnsemReal> x(dim);

    for(int j=0; j < dim; ++j)
    {
      x[j].resize(nbins);
      if (i == j)
	x[j] = Real(1);
      else
	x[j] = Real(0);
    }

    // Solve for the unit vector rhs
    Array<EnsemReal> col_vec = ensemSolveNxN(B, x);

    // Write onto each column
    for(int j=0; j < dim; ++j)
    {
      Binv[i][j] = col_vec[j];
    }
  }

  // Read REAL PART of eigenvectors
  Array2d<EnsemReal> A(dim,dim);
  for(int i=0; i < dim; ++i)
  {
    for(int j=0; j < dim; ++j)
    {
      ostringstream file_name;
      file_name << "ev_" << i << "_" << j << "_re.jknf";

      EnsemVectorReal tmp;
      read(file_name.str(), tmp);
      A[i][j] = peekObs(tmp, t);
    }
  }


  // form the IMAG PART of the inverse : b = inv( bb ) = - inv( A.Binv.A + B )
  Array2d<EnsemReal> bb(dim,dim);
  for(int i=0; i < dim; ++i)
  {
    for(int j=0; j < dim; ++j)
      {

	for(int k=0; k < dim; ++k)
	  { 
	    for(int l=0; l < dim; ++l)
	      {
		bb[i][j] +=  A[i][k]*Binv[k][l]*A[l][j];
	      }
	  }
	
	bb[i][j] += B[i][j];
	bb[i][j] = Real(-1)*bb[i][j];
      }
  }
  
  Array2d<EnsemReal> b(dim,dim);
  for(int i=0; i < dim; ++i)
    {
      Array<EnsemReal> xx(dim);
      
    for(int j=0; j < dim; ++j)
      {
      xx[j].resize(nbins);
      if (i == j)
	xx[j] = Real(1);
      else
	xx[j] = Real(0);
      }
    
    // Solve for the unit vector rhs
    Array<EnsemReal> col_vec = ensemSolveNxN(b, xx);

    // Write onto each column
    for(int j=0; j < dim; ++j)
      {
	b[i][j] = col_vec[j];
      }
    }

  // form the REAL PART of the inverse : a = -b * A * Binv
  Array2d<EnsemReal> a(dim,dim);
  for(int i=0; i < dim; ++i)
    {
      for(int j=0; j < dim; ++j)
	{
	  
	for(int k=0; k < dim; ++k)
	  { 
	    for(int l=0; l < dim; ++l)
	      {
		a[i][j] +=  b[i][k]*A[k][l]*Binv[l][j];
	      }
	  }
	
	a[i][j] = Real(-1)*a[i][j];
      }
    }

  // form the complex output matrix
  Array2d<EnsemComplex> Z(dim,dim);
  for(int i=0; i < dim; ++i)
  {
    for(int j=0; j < dim; ++j)
    {
      Z[i][j] = cmplx( a[i][j] , b[i][j] );
    }
  }

  // Write the inverse
  for(int i=0; i < dim; ++i)
  {
    for(int j=0; j < dim; ++j)
    {
      ostringstream file_name;
      file_name << "Z_" << i << "_" << j << "_t" << t <<".jknf";

      write(file_name.str(), Z[i][j]);
    }
  }
 

  return 0;
}
