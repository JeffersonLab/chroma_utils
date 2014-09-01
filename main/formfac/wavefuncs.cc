// $Id: wavefuncs.cc,v 2.0 2008/12/05 04:43:49 edwards Exp $
/*! \file
 * \brief Wavefunction overlaps
 */

#include "wavefuncs.h"
#include "formfac/formfac_qsq.h"
#include <iostream>

#include <math.h>

#define S3J_0		1e-10

#define S3J_MAX_FACT	25
#define S3J_TEST

#define S3J_EQUAL(a,b)		(fabs((a)-(b))<S3J_0)
#define S3J_MAX(a,b,c,ris)	(((a)>(b)?(ris=(a)):(ris=(b)))>(c)?ris:(ris=(c)))
#define S3J_MIN(a,b,c,ris)	(((a)<(b)?(ris=(a)):(ris=(b)))<(c)?ris:(ris=(c)))


using namespace std;

namespace FF
{
  int max(int a, int b) {return (a < b) ? b : a;}
  int min(int a, int b) {return (a < b) ? a : b;}


  //! The r-th unit vector
  /*! The r is 0-based */
  Array<Real> unitVec(int r)
  {
    Array<Real> unit(Nd-1);
    unit = zero;

    //    Real unit_norm = Real(1.0);

    if (r < 0 || r >= Nd-1)
    {
      cerr << __func__ << ": r is out of bounds" << endl;
      exit(1);
    }

    unit[r] = Real(1.0);

    return unit;
  }

  //! The r-th spherical basis vector
  /*! The r is 0-based */
  Array<Complex> spherVec(int r)
  {
    Array<Complex> vec(Nd-1);
    vec = zero;

    if (r < -1 || r > 1)
    {
      cerr << __func__ << ": r is out of bounds" << endl;
      exit(1);
    }
    else if(r == 0)
    {vec[2] = Real(1.0);}
    else if(r == 1)
    {vec[0] =  Real(-1.0 / sqrt(2.0) ); vec[1] = timesI( Real(-1.0 / sqrt(2.0)) ); }
    else if(r == -1)
    {vec[0] =  Real( 1.0 / sqrt(2.0) ); vec[1] = timesI( Real(-1.0 / sqrt(2.0)) ); }

    return vec;
  }


  //Wigner 3-j symbol
  double s3j(double j1, double j2, double j3, 
	     double m1, double m2, double m3) {
    
    /*  ( j1 j2 j3 )
	(            ) = delta(m1+m2+m3,0) * (-1)^(j1-j2-m3) * 
	( m1 m2 m3 )
	
	+-
	|  (j1+j2-j3)! (j1-j2+j3)! (-j1+j2+j3)! 
	* | -------------------------------------- ...
	|
	+-
	-+ 1/2
	(j1-m1)! (j1+m1)! (j2-m2)! (j2+m2)! (j3-m3)! (j3+m3)!  |
	... ------------------------------------------------------- |     * 
	(j1+j2+j3+1)!                       |
	-+
			      
	+---
	\                       (-1)^k
	*    |   ---------------------------------------------------------------------
	/      k! (j1+j2-j3-k)! (j1-m1-k)! (j2+m2-k)! (j3-j2+m1+k)! (j3-j1-m2+k)!
	+---
	k
			      
	Where factorials must have non-negative integral values:
			      
	j1+j2-j3   >= 0		j1-j2+j3   >= 0		-j1+j2+j3 >= 0		j1+j2+j3+1 >= 0
	k          >= 0		j1+j2-j3-k >= 0		j1-m1-k   >= 0		j2+m2-k    >= 0
	j3-j2+m1+k >= 0		j3-j1-m2+k >= 0
			      
	The 3j symbol is therefore non-null if
			      
	j1+j2 >= j3		(1)
	j1+j3 >= j2		(2)
	j2+j3 >= j1		(3)

	and k values in the sum must be such that
			      
	k <= j1+j2-j3		(4)			k >= 0				(7)
	k <= j1-m1			(5)			k >= -j3+j2-m1		(8)
	k <= j2+m2			(6)			k >= -j3+j1+m2		(9)
			      
	If no values of k satisfy the (4) to (9), the result is null because the sum is null,
	otherwise one can find   kmin < kmax   such that
  
	kmin <= k <= kmax
							
	(4) to (6) => kmin=MAX(j1+j2-j3,  j1-m1,       j2+m2     )
	(7) to (9) => kmax=MIN(0,         -j3+j2-m1,   -j3+j1+m2 )

	The condition kmin < kmax includes (1) to (3) because
							
	(4) and (7)    =>    (1)
	(5) and (8)    =>    (2)
	(6) and (9)    =>    (3)
		
	Once the values of kmin and kmax are found, the only "selection rule" is kmin<kmax.
    */
    
    int k, kmin, kmax;
    int jpm1, jmm1, jpm2, jmm2, jpm3, jmm3;
    int j1pj2mj3, j3mj2pm1, j3mj1mm2;
    double ris, mult, f[S3J_MAX_FACT];
    
    f[0]=1.0;
    mult=1.0;
    for (k=1; k<S3J_MAX_FACT; ++k) {
      
      f[k]=f[k-1]*mult;
      mult+=1.0;
    }
    
    jpm1=(int)(j1+m1);
    if (!S3J_EQUAL(jpm1,j1+m1)) return 0.0;
    jpm2=(int)(j2+m2);
    if (!S3J_EQUAL(jpm2,j2+m2)) return 0.0;
    jpm3=(int)(j3+m3);
    if (!S3J_EQUAL(jpm3,j3+m3)) return 0.0;
    jmm1=(int)(j1-m1);
    if (!S3J_EQUAL(jmm1,j1-m1)) return 0.0;
    jmm2=(int)(j2-m2);
    if (!S3J_EQUAL(jmm2,j2-m2)) return 0.0;
    jmm3=(int)(j3-m3);
    if (!S3J_EQUAL(jmm3,j3-m3)) return 0.0;
    
    /* delta(m1+m2+m3,0) */
    if ((jpm1-jmm1+jpm2-jmm2+jpm3-jmm3)!=0) return 0.0;
    
    /* j1+j2-j3 = (j1+j2-j3) + (m1+m2+m3) = jpm1+jpm2-jmm3 */
    j1pj2mj3=jpm1+jpm2-jmm3;
    /* j3-j2+m1 = (j3-j2+m1) - (m1+m2+m3) = jmm3-jpm2 */
    j3mj2pm1=jmm3-jpm2;
    /* j3-j1-m2 = (j3-j1-m2) + (m1+m2+m3) = jpm3-jmm1 */
    j3mj1mm2=jpm3-jmm1;

    S3J_MAX(-j3mj2pm1, -j3mj1mm2,   0,         kmin);
    S3J_MIN(j1pj2mj3,  jmm1,       jpm2,       kmax);
    if (kmin>kmax) return 0.0;	
    
    ris=0.0;
    if (kmin%2==0) mult=1.0;
    else mult=-1.0;
    for (k=kmin; k<=kmax; ++k) {
      ris+=mult/(f[k]*f[j1pj2mj3-k]*f[jmm1-k]*f[jpm2-k]*f[j3mj2pm1+k]*f[j3mj1mm2+k]);
      mult=-mult;
    }
    /* (-1)^(j1-j2-m3)=(-1)^(j1-j2-m3+m1+m2+m3)=(-1)^(jpm1-jmm2) */
    if ((jpm1-jmm2)%2!=0) ris=-ris;
    ris*=sqrt(f[j1pj2mj3]*f[jpm1-jmm2+jpm3]*f[-jmm1+jpm2+jpm3]*
	      f[jpm1]*f[jpm2]*f[jpm3]*f[jmm1]*f[jmm2]*f[jmm3]/
	      f[jpm1+jpm2+jpm3+1]);
    
    return ris;
  }
  
  //clebsch gordan coeffecient <j1 m1; j2 m2 | j3 m3>
  Real cleb(double j1, double m1, double j2, double m2, double j3, double m3)
  {
    double min_m3 = (-1) * m3;
    double phase = 1.0;
    int jjm = (int)(j1 - j2 + m3);
    if(jjm %2 != 0)
    {//odd 
      phase = -1.0;
    }
    double cleb_out = phase * sqrt( 2.0*j3 + 1.0 ) * s3j(j1, j2, j3, m1, m2, min_m3);
    return Real(cleb_out);
  }


  //! Energy from the continuum dispersion relation
  Real dispRel(const Real& mass, const ArrayInt& p, LatticeParam lattice)
  {
    int Ls = lattice.latt_size[0];
    const double xi  = lattice.aniso.xi;
    const double csq = lattice.aniso.c_sq;
    return sqrt(    ( Real(csq) / Real(xi*xi) ) * Real(norm2(contMom(p,Ls)))       + mass*mass);
  }

  // Make a Minkowski four vector momentum
  /* \return \f$p^{\mu}\f$*/
  Array<Real> make4Vec(const Real& mass, const ArrayInt& p, LatticeParam lattice)
  {
    Array<Real> vec(Nd);
    int Ls = lattice.latt_size[0];
    const double xi  = lattice.aniso.xi;
    const double csq = lattice.aniso.c_sq;

    vec[0] = Real( dispRel(mass, p, lattice) );

    Array<double> nnf = sqrt(csq) * contMom(p, Ls) /  xi;
    for(int i=0; i < 3; ++i)
      vec[i+1] = Real(nnf[i]);

    return vec;
  }


  //! Metric tensor
  /*! \f$ g^{\mu\nu} = g_{\mu\nu} = \{+1,-1,-1,-1\} \f$*/
  int metric(int mu, int nu)
  {
    if (mu < 0 || mu >= Nd)
    {
      cerr << __func__ << ": mu out of bounds" << endl;
      exit(1);
    }

    if (nu < 0 || nu >= Nd)
    {
      cerr << __func__ << ": nu out of bounds" << endl;
      exit(1);
    }

    if (mu != nu)
      return 0;

    // Now mu == nu
    if (mu == 0)
      return 1;
    else
      return -1;
  }



  //! 3D dot product
  Real dot3(const Array<Real>& p1, const Array<Real>& p2)
  {
    if ((p1.size() != p2.size()) && (p1.size() != Nd-1))
    {
      cerr << __func__ << ": not 4 vectors" << endl;
      exit(1);
    }

    Real sum = zero;

    for(int i=0; i < p1.size(); ++i)
      sum += p1[i] * p2[i];

    return sum;
  }

  //! 3D dot product
  Complex dot3(const Array<Complex>& p1, const Array<Real>& p2)
  {
    if ((p1.size() != p2.size()) && (p1.size() != Nd-1))
    {
      cerr << __func__ << ": not 4 vectors" << endl;
      exit(1);
    }

    Complex sum = zero;

    for(int i=0; i < p1.size(); ++i)
      sum += p1[i] * p2[i];

    return sum;
  }

  //! Minkowski 4D dot product
  Real dot4(const Array<Real>& p1, const Array<Real>& p2)
  {
    if ((p1.size() != p2.size()) && (p1.size() != Nd))
    {
      cerr << __func__ << ": not 4 vectors" << endl;
      exit(1);
    }

    Real sum = p1[0] * p2[0];

    for(int i=1; i < p1.size(); ++i)
      sum -= p1[i] * p2[i];

    return sum;
  }


  //! Raise 4-vector
  /*! \return \f$g^{\mu\nu}*p_{\nu}\f$*/
  template<typename T>
  Array<T> raise4Vec(const Array<T>& p)
  {
    if (p.size() != Nd)
    {
      cerr << __func__ << ": not a 4 vector" << endl;
      exit(1);
    }

    Array<T> vec(Nd);

    for(int mu=0; mu < vec.size(); ++mu)
      vec[mu] = Real(metric(mu,mu))*p[mu];

    return vec;
  }


  //! Lower 4-vector
  /*! \return \f$g_{\mu\nu}*p^{\nu}\f$*/
  template<typename T>
  Array<T> lower4Vec(const Array<T>& p)
  {
    return raise4Vec(p);
  }


#if 0
  //! Minkowski polarization vector
  /*! \return \f$\epsilon^{\mu}\f$*/
  Array<Real> minkPolVec(const Real& mass, const ArrayInt& p, int r, LatticeParam lattice)
  {
    Array<Real> unit = unitVec(r);

    int Ls = lattice.latt_size[0];
    const double xi  = lattice.aniso.xi;
    const double csq = lattice.aniso.c_sq;

    Array<double> nnf = sqrt(csq) * contMom(p, Ls) / xi;
    Array<Real> nf(nnf.size());
    for(int i=0; i < nf.size(); ++i)
      nf[i] = nnf[i];

    Array<Real> pol(Nd);
    Real dott = dot3(unit, nf) / mass;
    pol[0] = dott;

    for(int k=0; k < Nd-1; ++k)
      pol[k+1] = unit[k] + dott * (nf[k] / (dispRel(mass,p, lattice) + mass));

    return pol;
  }
#endif

  //! Minkowski polarization vector
  /*! \return \f$\epsilon^{\mu}\f$*/
  Array<Complex> minkPolVec(const Real& mass, const ArrayInt& p, int r, LatticeParam lattice)
  {
    Array<Complex> unit = spherVec(r);

    int Ls = lattice.latt_size[0];
    const double xi  = lattice.aniso.xi;
    const double csq = lattice.aniso.c_sq;

    Array<double> nnf = sqrt(csq) * contMom(p, Ls) / xi;
    Array<Real> nf(nnf.size());
    for(int i=0; i < nf.size(); ++i)
      nf[i] = nnf[i];

    Array<Complex> pol(Nd);
    Complex dott = dot3(unit, nf) / mass;
    pol[0] = dott;

    for(int k=0; k < Nd-1; ++k)
      pol[k+1] = unit[k] + dott * (nf[k] / (dispRel(mass,p, lattice) + mass));

    return pol;
  }

#if 0
  //! Minkowski rank-2 polarization tensor
  /*! 
   *
   * \return \f$\epsilon^{\mu\nu}\f$
   */
  Array< Array<Real> > minkPolTens(const Real& mass, const ArrayInt& p, int r, LatticeParam lattice)
  {

    Array< Array<Real> >  tens(Nd);
    for(int mu=0; mu < Nd; ++mu)
      tens[mu].resize(Nd);
   
    Array<Real> eps0 = minkPolVec(mass,p,0,lattice);
    Array<Real> eps1 = minkPolVec(mass,p,1,lattice);
    Array<Real> eps2 = minkPolVec(mass,p,2,lattice);

    for(int mu=0; mu < Nd; ++mu)
    {
      for(int nu=0; nu < Nd; ++nu)
      {
	if(r==0)
	{tens[mu][nu] = (eps1[mu]*eps0[nu] + eps0[mu]*eps1[nu]) /sqrt(Real(2.0));   }	    
	else if(r==1)
	{tens[mu][nu] = (eps2[mu]*eps0[nu] + eps0[mu]*eps2[nu]) /sqrt(Real(2.0));   }
	else if(r==2)
	{tens[mu][nu] = (eps2[mu]*eps1[nu] + eps1[mu]*eps2[nu])  /sqrt(Real(2.0));   }
	else if(r==3)
	{tens[mu][nu] = (eps0[mu]*eps0[nu] + eps1[mu]*eps1[nu]  - Real(2.0)* eps2[mu]*eps2[nu] ) /sqrt(Real(6.0));   }
	else if(r==4)
	{tens[mu][nu] = ( eps0[mu]*eps0[nu] - eps1[mu]*eps1[nu] ) /sqrt(Real(2.0));   };
      }
    }

    return tens;
  }
#endif

  //! Minkowski rank-2 polarization tensor
  /*! 
   *
   * \return \f$\epsilon^{\mu\nu}\f$
   */
  Array< Array<Complex> > minkPolTens(const Real& mass, const ArrayInt& p, int r, LatticeParam lattice)
  {

    Array< Array<Complex> >  tens(Nd);
    for(int mu=0; mu < Nd; ++mu)
      tens[mu].resize(Nd);

    for (int mu=0; mu < Nd; ++mu)
      for (int nu=0; nu < Nd; ++nu)
      {tens[mu][nu] = zero;}
   
    Array< Array<Complex> > eps(3);
    for (int rr=0; rr < 3; ++rr)
      eps[rr].resize(Nd);
    //eps[r][mu]
    for (int rr=0; rr < 3; ++rr)
    {	eps[rr] = minkPolVec(mass, p, rr - 1 , lattice);}

    for (int mu=0; mu < Nd; ++mu)
      for (int nu=0; nu < Nd; ++nu)
	for (int rr1 = -1; rr1 < 2; ++rr1)
	  for (int rr2 = -1; rr2 < 2; ++rr2)
	  { tens[mu][nu] += cleb(1.0, toFloat( rr1 ), 1.0, toFloat( rr2 ), 2.0, toFloat( r ) ) * eps[rr1 + 1][mu] * eps[rr2 + 1][nu];}

    return tens;
  }

  //! Minkowski rank-3 polarization tensor
  /*! 
   *
   * \return \f$\epsilon^{\mu\nu\tau}\f$
   */
  Array< Array< Array<Complex> > > minkPolTensSpin3(const Real& mass, const ArrayInt& p, int r, LatticeParam lattice)
  {
    Array< Array< Array<Complex> > >  tens(Nd);
    for(int mu=0; mu < Nd; ++mu)
      tens[mu].resize(Nd);
    for(int mu=0; mu < Nd; ++mu)
      for(int nu=0; nu < Nd; ++nu)
      {
	tens[mu][nu].resize(Nd);
      }
    for(int mu=0; mu < Nd; ++mu)
      for(int nu=0; nu < Nd; ++nu)
	for(int tau=0; tau < Nd; ++tau)
	  tens[mu][nu][tau] = zero;

    Array< Array<Complex> > eps(3);
    for (int rr=0; rr < 3; ++rr)
      eps[rr].resize(Nd);
    //eps[r][mu]
    for (int rr=0; rr < 3; ++rr)
    {	eps[rr] = minkPolVec(mass, p, rr - 1 , lattice);}

    Array< Array< Array<Complex> > > eps2(5);
    for (int rr=0; rr < 5; ++rr)
    {
      eps2[rr].resize(Nd);
      for (int mu=0; mu< Nd; ++mu)
	eps2[rr][mu].resize(Nd);
    }
    //eps2[r][mu][nu]
    for (int rr=0; rr < 5; ++rr)
    {	eps2[rr] = minkPolTens(mass, p, rr - 2 , lattice);}

    for (int mu=0; mu < Nd; ++mu)
      for (int nu=0; nu < Nd; ++nu)
	for (int tau=0; tau < Nd; ++tau)
	  for (int rr1 = -1; rr1 < 2; ++rr1)
	    for (int rr2 = -2; rr2 < 3; ++rr2)
	    { tens[mu][nu][tau] += cleb(1.0, toFloat( rr1 ), 2.0, toFloat( rr2 ), 3.0, toFloat( r ) ) * eps[rr1 + 1][mu] * eps2[rr2 + 2][nu][tau];}
   
    return tens;
  }


  //---------------------------------------------------------------------------

  template<typename T> //so this could be an array or anything else
  class Wvf0index
  {
  public:
    virtual Array<T> operator()(void) const = 0;
    virtual int getNumTerms() const = 0;
  };

  class WvfA0 : public Wvf0index<Complex>
  {
  public:
    //constructor 
    WvfA0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_), lattice(lattice_) {}
    
    Array<Complex> operator()(void) const
      {
	Array<Complex> out(1);
	out[0] = Real(1);
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };

  class WvfGamma5Pi0 : public Wvf0index<Complex>
  {
  public:
    //constructor 
    WvfGamma5Pi0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_), lattice(lattice_) {}
    
    Array<Complex> operator()(void) const
      {
	Array<Complex> out(1);
	out[0] = timesI( Real(1) );
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };

  //--------------------------------------------------------------------------
  Array<Complex> test_wvf0index(const Real& mass, const ArrayInt& p, LatticeParam lattice)
  {
    WvfA0 wvf(mass, p, lattice);
    return wvf();
  }


  //-------------------------------------------------------------------------------------
  template<typename T> //so this could be an array
  class Wvf1index
  {
  public:
    virtual Array<T> operator()(int mu) const = 0;
    virtual int getNumTerms() const = 0;
  };

  class WvfGammaMuB0 : public Wvf1index<Complex>
  {
  public:
    //constructor 
    WvfGammaMuB0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_), lattice(lattice_) {}
    
    Array<Complex> operator()(int mu) const
      {
	Array<Real> pp = make4Vec(mass,p,lattice);
	Array<Complex> out(1);
	out[0] = timesI(pp[mu]);
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };

  class WvfGammaMuRho1 : public Wvf1index<Complex>
  {
  public:
    //constructor 
    WvfGammaMuRho1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_), lattice(lattice_) {}
    
    Array<Complex> operator()(int mu) const
      {
	Array<Complex> eps = minkPolVec(mass,p,r,lattice);
	Array<Complex> out(1);
	out[0] = eps[mu];
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfGamma5GammaMuPi0 : public Wvf1index<Complex>
  {
  public:
    //constructor 
    WvfGamma5GammaMuPi0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu) const
      {
	Array<Real> pp = make4Vec(mass,p,lattice);
	Array<Complex> out(1);
	out[0] = timesI(pp[mu]);
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };

  class WvfGamma5GammaMuA1 : public Wvf1index<Complex>
  {
  public:
    //constructor 
    WvfGamma5GammaMuA1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu) const
      {
	Array<Complex> eps = minkPolVec(mass,p,r,lattice);
	Array<Complex> out(1);
	out[0] = eps[mu];
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfNablaMuB0 : public Wvf1index<Complex>
  {
  public:
    //constructor 
    WvfNablaMuB0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_), lattice(lattice_) {}
    
    Array<Complex> operator()(int mu) const
      {
	Array<Real> pp = make4Vec(mass,p,lattice);
	Array<Complex> out(1);
	out[0] = pp[mu];
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };

  class WvfNablaMuRho1 : public Wvf1index<Complex>
  {
  public:
    //constructor 
    WvfNablaMuRho1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_), lattice(lattice_) {}
    
    Array<Complex> operator()(int mu) const
      {
	Array<Complex> eps = minkPolVec(mass,p,r,lattice);
	Array<Complex> out(1);
	out[0] = timesI(eps[mu]);
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfGamma5NablaMuRho0 : public Wvf1index<Complex>
  {
  public:
    //constructor 
    WvfGamma5NablaMuRho0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu) const
      {
	Array<Real> pp = make4Vec(mass,p,lattice);
	Array<Complex> out(1);
	out[0] = timesI(pp[mu]);
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };

  class WvfGamma5NablaMuB1 : public Wvf1index<Complex>
  {
  public:
    //constructor 
    WvfGamma5NablaMuB1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu) const
      {
	Array<Complex> eps = minkPolVec(mass,p,r,lattice);
	Array<Complex> out(1);
	out[0] = timesI(eps[mu]);
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  //--------------------------------------------------------------------------
  Array<Complex> test_wvf1index(const Real& mass, const ArrayInt& p, int r, LatticeParam lattice, int mu)
  {
    WvfGammaMuRho1 wvf(mass, p, r, lattice);
    return wvf(mu);
  }



  //----------------------------------------------------------------------------------
  template<typename T> //so this could be an array
  class Wvf2index
  {
  public:
    virtual Array<T> operator()(int mu, int nu) const = 0;
    virtual int getNumTerms() const = 0;
  };



  class WvfSigmaMuNuB1 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfSigmaMuNuB1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array<Real> pp = lower4Vec(make4Vec(mass,p,lattice));

	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 

	Array<Complex> eps = lower4Vec(minkPolVec(mass,p,r,lattice));
	Complex sum = zero;
	
	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int delta=0; delta < Nd; ++delta)
	  {
	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = delta;
	      
	    Real a = peekObs(anti,n1);
	    sum += a * eps[alpha] * pp[delta];
	  }
	Array<Complex> out(1);
	out[0] =  sum; 
	return out;
      };

    int getNumTerms() const
      {
	return 1;
      }


  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  Real test_antisym()
  {
    TensorReal anti = antiSymTensor<4>();
    Array<int> n1(4);
    n1[0] = 0; n1[1] = 1; n1[2] = 2; n1[3] = 3;
    Real a = peekObs(anti,n1);
    return a;
  }


  class WvfSigmaMuNuRho1 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfSigmaMuNuRho1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array<Real> pp = make4Vec(mass,p,lattice);
	Array<Complex> eps = minkPolVec(mass,p,r,lattice);

	Array<Complex> out(1);
	out[0] = timesI(eps[mu]*pp[nu] - eps[nu]*pp[mu]); 
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }

  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfGammaMuNablaNuA0 : public Wvf2index<Complex >
  {
  public:
    //constructor 
    WvfGammaMuNablaNuA0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array<Complex> out(2);
	out[0] = timesI( Real(metric(mu,nu)) );
	Array<Real> pp = make4Vec(mass,p,lattice);
	out[1] =  timesI( pp[mu]*pp[nu]  );
	return out;
      };
    int getNumTerms() const
      {
	return 2;
      }

  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };


  class WvfGammaMuNablaNuA1 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfGammaMuNablaNuA1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
        Array<Real> pp = lower4Vec(make4Vec(mass,p,lattice));

	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 

	Array<Complex> eps = lower4Vec(minkPolVec(mass,p,r,lattice));
	Complex sum = zero;

	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int delta=0; delta < Nd; ++delta)
	  {
	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = delta;
	
	    Real a = peekObs(anti,n1);
	    sum += a * eps[alpha] * pp[delta];
	  };
	Array<Complex> out(1);
	out[0] = sum;
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };


  class WvfGammaMuNablaNuA2 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfGammaMuNablaNuA2(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array< Array<Complex> > tens = minkPolTens(mass,p,r,lattice);
	Array<Complex> out(1);
	out[0] =  timesI(tens[mu][nu]);
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  //--------------------------------------------------------------------------
  Array<Complex> test_wvf2index(const Real& mass, const ArrayInt& p, int r, LatticeParam lattice, int mu, int nu)
  {
    WvfGammaMuNablaNuA2 wvf(mass, p,r, lattice);
    return wvf(mu, nu);
  }
  //-------------------------------------------------------------------------------

  class WvfGammaMuNablaNuPi1 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfGammaMuNablaNuPi1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array<Real> pp = make4Vec(mass,p,lattice);
	Array<Complex> eps = minkPolVec(mass,p,r,lattice);
	Array<Complex> out(2);
	out[0] =  timesI(eps[mu]*pp[nu] + eps[nu]*pp[mu]);
	out[1] =  timesI(eps[mu]*pp[nu] - eps[nu]*pp[mu]);
	return out;
      };
    int getNumTerms() const
      {
	return 2;
      }

  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };


  class WvfGamma5GammaMuNablaNuRho0 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfGamma5GammaMuNablaNuRho0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array<Complex> out(2);
	out[0] = timesI( Real(metric(mu,nu)) );
	Array<Real> pp = make4Vec(mass,p,lattice);
	out[1] =  timesI( pp[mu]*pp[nu]  );
	return out;
      };
    int getNumTerms() const
      {
	return 2;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };


  class WvfGamma5GammaMuNablaNuRho1 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfGamma5GammaMuNablaNuRho1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
        Array<Real> pp = lower4Vec(make4Vec(mass,p,lattice));

	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 

	Array<Complex> eps = lower4Vec(minkPolVec(mass,p,r,lattice));
	Complex sum = zero;

	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int delta=0; delta < Nd; ++delta)
	  {
	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = delta;
	
	    Real a = peekObs(anti,n1);
	    sum += a * eps[alpha] * pp[delta];
	  }
	Array<Complex> out(1);
	out[0] =  sum;
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };


  class WvfGamma5GammaMuNablaNuRho2 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfGamma5GammaMuNablaNuRho2(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array< Array<Complex> > tens = minkPolTens(mass,p,r,lattice);
	Array<Complex> out(1);
	out[0] =  timesI(tens[mu][nu]);
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };


  class WvfGamma5GammaMuNablaNuB1 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfGamma5GammaMuNablaNuB1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array<Real> pp = make4Vec(mass,p,lattice);
	Array<Complex> eps = minkPolVec(mass,p,r,lattice);
	Array<Complex> out(2);
	out[0] =  timesI(eps[mu]*pp[nu] + eps[nu]*pp[mu]);
	out[1] =  timesI(eps[mu]*pp[nu] - eps[nu]*pp[mu]);
	return out;
      };
    int getNumTerms() const
      {
	return 2;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };
 
  class WvfNablaMuNablaNuA0 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfNablaMuNablaNuA0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array<Complex> out(2);
	out[0] = Real(metric(mu,nu)) ;
	Array<Real> pp = make4Vec(mass,p,lattice);
	out[1] =   pp[mu]*pp[nu] ;
	return out;
      };
    int getNumTerms() const
      {
	return 2;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };



  class WvfNablaMuNablaNuA1 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfNablaMuNablaNuA1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
        Array<Real> pp = lower4Vec(make4Vec(mass,p,lattice));

	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 

	Array<Complex> eps = lower4Vec(minkPolVec(mass,p,r,lattice));
	Complex sum = zero;

	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int delta=0; delta < Nd; ++delta)
	  {
	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = delta;
	
	    Real a = peekObs(anti,n1);
	    sum += a * eps[alpha] * pp[delta];
	  }
	Array<Complex> out(1);
	out[0] = timesI(sum); 
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfNablaMuNablaNuA2 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfNablaMuNablaNuA2(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array< Array<Complex> > tens = minkPolTens(mass,p,r,lattice);
	Array< Complex> out(1);
	out[0] =  tens[mu][nu]; 
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };


  class WvfNablaMuNablaNuPi1 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfNablaMuNablaNuPi1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array<Real> pp = make4Vec(mass,p,lattice);
	Array<Complex> eps = minkPolVec(mass,p,r,lattice);
	Array<Complex> out(2);
	out[0] =  eps[mu]*pp[nu] + eps[nu]*pp[mu];
	out[1] =  eps[mu]*pp[nu] - eps[nu]*pp[mu];
	return out;
      };
    int getNumTerms() const
      {
	return 2;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfGamma5NablaMuNablaNuPi0 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfGamma5NablaMuNablaNuPi0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array<Complex> out(2);
	out[0] = timesI( Real(metric(mu,nu)) );
	Array<Real> pp = make4Vec(mass,p,lattice);
	out[1] =  timesI( pp[mu]*pp[nu]  );
	return out;
      };
    int getNumTerms() const
      {
	return 2;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };


  class WvfGamma5NablaMuNablaNuPi1 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfGamma5NablaMuNablaNuPi1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
        Array<Real> pp = lower4Vec(make4Vec(mass,p,lattice));

	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 

	Array<Complex> eps = lower4Vec(minkPolVec(mass,p,r,lattice));
	Complex sum = zero;

	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int delta=0; delta < Nd; ++delta)
	  {
	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = delta;
	
	    Real a = peekObs(anti,n1);
	    sum += a * eps[alpha] * pp[delta];
	  }
	Array<Complex> out(1);
	out[0] =  timesI(sum); 
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };


  class WvfGamma5NablaMuNablaNuPi2 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfGamma5NablaMuNablaNuPi2(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array< Array<Complex> > tens = minkPolTens(mass,p,r,lattice);
	Array<Complex> out(1);
	out[0] = tens[mu][nu]; 
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };


  class WvfGamma5NablaMuNablaNuA1 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfGamma5NablaMuNablaNuA1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array<Real> pp = make4Vec(mass,p,lattice);
	Array<Complex> eps = minkPolVec(mass,p,r,lattice);
	Array<Complex> out(2);
	out[0] =  eps[mu]*pp[nu] + eps[nu]*pp[mu];
	out[1] =  eps[mu]*pp[nu] - eps[nu]*pp[mu];
	return out;
      };
    int getNumTerms() const
      {
	return 2;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };


  class WvfFMuNuB0 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfFMuNuB0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array<Complex> out(2);
	out[0] = timesI( Real(metric(mu,nu)) );
	Array<Real> pp = make4Vec(mass,p,lattice);
	out[1] =  timesI( pp[mu]*pp[nu]  );
	return out;
      };
    int getNumTerms() const
      {
	return 2;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };


  class WvfFMuNuB1 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfFMuNuB1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
        Array<Real> pp = lower4Vec(make4Vec(mass,p,lattice));

	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 

	Array<Complex> eps = lower4Vec(minkPolVec(mass,p,r,lattice));
	Complex sum = zero;

	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int delta=0; delta < Nd; ++delta)
	  {
	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = delta;
	
	    Real a = peekObs(anti,n1);
	    sum += a * eps[alpha] * pp[delta];
	  }

	Array<Complex> out(1);
	out[0] = timesI(sum); 
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };


  class WvfFMuNuB2 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfFMuNuB2(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array< Array<Complex> > tens = minkPolTens(mass,p,r,lattice);
	Array<Complex> out(1);
	out[0] =  tens[mu][nu]; 
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };


  class WvfFMuNuRho1 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfFMuNuRho1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array<Real> pp = make4Vec(mass,p,lattice);
	Array<Complex> eps = minkPolVec(mass,p,r,lattice);
	Array<Complex> out(1);
	out[0] =  eps[mu]*pp[nu] - eps[nu]*pp[mu];
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };


  class WvfGamma5FMuNuRho0 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfGamma5FMuNuRho0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array<Complex> out(2);
	out[0] = Real(metric(mu,nu)) ;
	Array<Real> pp = make4Vec(mass,p,lattice);
	out[1] =  pp[mu]*pp[nu]  ;
	return out;
      };
    int getNumTerms() const
      {
	return 2;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };


  class WvfGamma5FMuNuRho1 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfGamma5FMuNuRho1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
        Array<Real> pp = lower4Vec(make4Vec(mass,p,lattice));

	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 

	Array<Complex> eps = lower4Vec(minkPolVec(mass,p,r,lattice));
	Complex sum = zero;

	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int delta=0; delta < Nd; ++delta)
	  {
	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = delta;
	
	    Real a = peekObs(anti,n1);
	    sum += a * eps[alpha] * pp[delta];
	  }
	Array<Complex> out(1);
	out[0]  =  timesI(sum);
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };


  class WvfGamma5FMuNuRho2 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfGamma5FMuNuRho2(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array< Array<Complex> > tens = minkPolTens(mass,p,r,lattice);
	Array<Complex> out(1);
	out[0] = tens[mu][nu]; 
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };


  class WvfGamma5FMuNuB1 : public Wvf2index<Complex>
  {
  public:
    //constructor 
    WvfGamma5FMuNuB1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu) const
      {
	Array<Real> pp = make4Vec(mass,p,lattice);
	Array<Complex> eps = minkPolVec(mass,p,r,lattice);
	Array<Complex> out(1);
	out[0] =  eps[mu]*pp[nu] - eps[nu]*pp[mu];
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };



  //----------------------------------------------------------------------------------
  template<typename T> //so this could be an array
  class Wvf3index
  {
  public:
    virtual Array<T> operator()(int mu, int nu, int tau) const = 0;
    virtual int getNumTerms() const = 0;
  };




  class WvfGammaMuGammaNuNablaTauA0 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGammaMuGammaNuNablaTauA0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Real> pp = make4Vec(mass,p,lattice);
	Array<Complex> out(4);  
	out[0] =  timesI( Real(metric(mu,tau)) * pp[nu] +  Real(metric(nu,tau)) * pp[mu] )   ;
	out[1] =  timesI( Real(metric(mu,tau)) * pp[nu] -  Real(metric(nu,tau)) * pp[mu] )   ;
	out[2] =  timesI( Real(metric(mu,nu)) * pp[tau] ) ;
	out[3] =  timesI( pp[mu] * pp[nu] * pp[tau] );

	return out;
      };
    int getNumTerms() const
      {
	return 4;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };


  class WvfGammaMuGammaNuNablaTauPi0 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGammaMuGammaNuNablaTauPi0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Real> pp = lower4Vec( make4Vec(mass,p,lattice) );
	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 
	Complex sum = zero;

	for(int alpha=0; alpha < Nd; ++alpha)
	{
	  n1[0] = mu; n1[1] = nu; n1[2] = tau; n1[3] = alpha;
	
	  Real a = peekObs(anti,n1);
	  sum += a * pp[alpha];
	}
	Array<Complex> out(1);
	out[0] = timesI(sum);
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };

  class WvfGammaMuGammaNuNablaTauA1 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGammaMuGammaNuNablaTauA1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(4);
	Array<Complex> eps = lower4Vec(minkPolVec(mass,p,r,lattice));
	Array<Real> pp = make4Vec(mass,p,lattice);
	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 
	for(int i=0; i < 4; ++i)	{out[i] = zero;}

	for(int alpha=0; alpha < Nd; ++alpha)
	{
	  n1[0] = mu; n1[1] = nu; n1[2] = tau; n1[3] = alpha;
	
	  Real a = peekObs(anti,n1);
	  out[0] += a * eps[alpha];
	}
	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int beta=0; beta < Nd; ++beta) 
	  {
	    n1[0] = mu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;
	    Real a = peekObs(anti,n1);
	    Complex d1 = a * eps[beta] * lower4Vec(pp)[alpha] * pp[nu];	    
	   
	    n1[0] = nu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;
	    a = peekObs(anti,n1);
	    Complex d2 = a * eps[beta] * lower4Vec(pp)[alpha] * pp[mu];	 	   
	    
	    out[1] += d1 + d2;
	    out[2] += d1 - d2;

	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = beta;
	    a = peekObs(anti,n1);
	    out[3] += a * eps[beta] * lower4Vec(pp)[alpha] * pp[tau];
	  }

	return out;
      };
    int getNumTerms() const
      {
	return 4;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfGammaMuGammaNuNablaTauPi1 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGammaMuGammaNuNablaTauPi1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(6);
	Array<Complex> eps = minkPolVec(mass,p,r,lattice);
	Array<Real> pp = make4Vec(mass,p,lattice);

	out[0] = timesI( Real(metric(mu,tau)) * eps[nu] +  Real(metric(nu,tau)) * eps[mu]);
	out[1] = timesI( Real(metric(mu,tau)) * eps[nu] -  Real(metric(nu,tau)) * eps[mu]);
	out[2] = timesI( Real(metric(mu,nu)) * eps[tau]);
	out[3] = timesI( pp[mu]* pp[tau] * eps[nu] + pp[nu] * pp[tau] * eps[mu]);
	out[4] = timesI( pp[mu]* pp[tau] * eps[nu] - pp[nu] * pp[tau] * eps[mu]);
	out[5] = timesI( pp[mu] * pp[nu] * eps[tau]);

	return out;
      };
    int getNumTerms() const
      {
	return 6;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfGammaMuGammaNuNablaTauA2 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGammaMuGammaNuNablaTauA2(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(3);
	Array< Array<Complex> > tens = minkPolTens(mass,p,r,lattice); 
	Array<Real> pp = make4Vec(mass,p,lattice);

	out[0] =  timesI( tens[mu][tau] * pp[nu] +  tens[nu][tau] * pp[mu] );
	out[1] =  timesI( tens[mu][tau] * pp[nu] -  tens[nu][tau] * pp[mu] );
	out[2] =  timesI( tens[mu][nu] * pp[tau] );

	return out;
      };
    int getNumTerms() const
      {
	return 3;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfGammaMuGammaNuNablaTauPi2 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGammaMuGammaNuNablaTauPi2(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(3);
	for(int i=0; i < 3; ++i)	{out[i] = zero;}

	Array< Array<Complex> > tens = minkPolTens(mass,p,r,lattice); 
	Array<Real> pp = lower4Vec( make4Vec(mass,p,lattice) );

	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 

	Array<Complex> tempmu = lower4Vec( tens[mu] );
	Array<Complex> tempnu = lower4Vec( tens[nu] );
	Array<Complex> temptau = lower4Vec( tens[tau] );

	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int beta=0; beta < Nd; ++beta) 
	  {
	    n1[0] = mu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;	    
	    Real a = peekObs(anti,n1);
	    Complex d1 = a * tempnu[alpha] * pp[beta] ;

	    n1[0] = nu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;	    
	    a = peekObs(anti,n1);
	    Complex d2 = a * tempmu[alpha] * pp[beta] ;
	    
	    out[0] += d1 + d2;
	    out[1] += d1 - d2;

	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = beta;	    
	    a = peekObs(anti,n1);
	    out[2] = a * temptau[alpha] * pp[beta] ;	    
	  }  

	return out;
      };
    int getNumTerms() const
      {
	return 3;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfGammaMuGammaNuNablaTauPi3 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGammaMuGammaNuNablaTauPi3(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array< Array< Array<Complex> > > tens3 = minkPolTensSpin3(mass,p,r,lattice); 
	Array<Complex> out(1);
	out[0] = timesI( tens3[mu][nu][tau] );

	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };


  // Nabla Nabla Gamma


  class WvfNablaMuNablaNuGammaTauB0 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfNablaMuNablaNuGammaTauB0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Real> pp = make4Vec(mass,p,lattice);
	Array<Complex> out(4);  
	out[0] =  timesI( Real(metric(mu,tau)) * pp[nu] +  Real(metric(nu,tau)) * pp[mu] )   ;
	out[1] =  timesI( Real(metric(mu,tau)) * pp[nu] -  Real(metric(nu,tau)) * pp[mu] )   ;
	out[2] =  timesI( Real(metric(mu,nu)) * pp[tau] ) ;
	out[3] =  timesI( pp[mu] * pp[nu] * pp[tau] );

	return out;
      };
    int getNumTerms() const
      {
	return 4;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };

  class WvfNablaMuNablaNuGammaTauRho0 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfNablaMuNablaNuGammaTauRho0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Real> pp = lower4Vec( make4Vec(mass,p,lattice) );
	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 
	Complex sum = zero;

	for(int alpha=0; alpha < Nd; ++alpha)
	{
	  n1[0] = mu; n1[1] = nu; n1[2] = tau; n1[3] = alpha;
	
	  Real a = peekObs(anti,n1);
	  sum += a * pp[alpha];
	}
	Array<Complex> out(1);
	out[0] = timesI(sum);
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };

  class WvfNablaMuNablaNuGammaTauB1 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfNablaMuNablaNuGammaTauB1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(4);
	Array<Complex> eps = lower4Vec(minkPolVec(mass,p,r,lattice));
	Array<Real> pp = make4Vec(mass,p,lattice);
	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 
	for(int i=0; i < 4; ++i)	{out[i] = zero;}

	for(int alpha=0; alpha < Nd; ++alpha)
	{
	  n1[0] = mu; n1[1] = nu; n1[2] = tau; n1[3] = alpha;
	
	  Real a = peekObs(anti,n1);
	  out[0] += timesI( a * eps[alpha] );
	}
	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int beta=0; beta < Nd; ++beta) 
	  {
	    n1[0] = mu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;
	    Real a = peekObs(anti,n1);
	    Complex d1 = a * eps[beta] * lower4Vec(pp)[alpha] * pp[nu];	    
	   
	    n1[0] = nu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;
	    a = peekObs(anti,n1);
	    Complex d2 = a * eps[beta] * lower4Vec(pp)[alpha] * pp[mu];	 	   
	    
	    out[1] += timesI(d1 + d2);
	    out[2] += timesI(d1 - d2);

	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = beta;
	    a = peekObs(anti,n1);
	    out[3] += timesI( a * eps[beta] * lower4Vec(pp)[alpha] * pp[tau] );
	  }

	return out;
      };
    int getNumTerms() const
      {
	return 4;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfNablaMuNablaNuGammaTauRho1 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfNablaMuNablaNuGammaTauRho1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(6);
	Array<Complex> eps = minkPolVec(mass,p,r,lattice);
	Array<Real> pp = make4Vec(mass,p,lattice);

	out[0] = Real(metric(mu,tau)) * eps[nu] +  Real(metric(nu,tau)) * eps[mu];
	out[1] = Real(metric(mu,tau)) * eps[nu] -  Real(metric(nu,tau)) * eps[mu];
	out[2] = Real(metric(mu,nu)) * eps[tau];
	out[3] = pp[mu]* pp[tau] * eps[nu] + pp[nu] * pp[tau] * eps[mu];
	out[4] = pp[mu]* pp[tau] * eps[nu] - pp[nu] * pp[tau] * eps[mu];
	out[5] = pp[mu] * pp[nu] * eps[tau];

	return out;
      };
    int getNumTerms() const
      {
	return 6;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfNablaMuNablaNuGammaTauB2 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfNablaMuNablaNuGammaTauB2(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(3);
	Array< Array<Complex> > tens = minkPolTens(mass,p,r,lattice); 
	Array<Real> pp = make4Vec(mass,p,lattice);

	out[0] = tens[mu][tau] * pp[nu] +  tens[nu][tau] * pp[mu];
	out[1] = tens[mu][tau] * pp[nu] -  tens[nu][tau] * pp[mu];
	out[2] = tens[mu][nu] * pp[tau];

	return out;
      };
    int getNumTerms() const
      {
	return 3;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfNablaMuNablaNuGammaTauRho2 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfNablaMuNablaNuGammaTauRho2(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(3);
	for(int i=0; i < 3; ++i)	{out[i] = zero;}

	Array< Array<Complex> > tens = minkPolTens(mass,p,r,lattice); 
	Array<Real> pp = lower4Vec( make4Vec(mass,p,lattice) );

	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 

	Array<Complex> tempmu = lower4Vec( tens[mu] );
	Array<Complex> tempnu = lower4Vec( tens[nu] );
	Array<Complex> temptau = lower4Vec( tens[tau] );

	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int beta=0; beta < Nd; ++beta) 
	  {
	    n1[0] = mu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;	    
	    Real a = peekObs(anti,n1);
	    Complex d1 = a * tempnu[alpha] * pp[beta] ;

	    n1[0] = nu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;	    
	    a = peekObs(anti,n1);
	    Complex d2 = a * tempmu[alpha] * pp[beta] ;
	    
	    out[0] += timesI(d1 + d2);
	    out[1] += timesI(d1 - d2);

	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = beta;	    
	    a = peekObs(anti,n1);
	    out[2] = timesI( a * temptau[alpha] * pp[beta] ) ;	    
	  }  

	return out;
      };
    int getNumTerms() const
      {
	return 3;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfNablaMuNablaNuGammaTauRho3 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfNablaMuNablaNuGammaTauRho3(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array< Array< Array<Complex> > > tens3 = minkPolTensSpin3(mass,p,r,lattice); 
	Array<Complex> out(1);
	out[0] = tens3[mu][nu][tau];

	return out;

      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };



  //Gamma5 NablaMu NablaNu GammaTau


  class WvfGamma5NablaMuNablaNuGammaTauPi0 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGamma5NablaMuNablaNuGammaTauPi0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Real> pp = make4Vec(mass,p,lattice);
	Array<Complex> out(4);  
	out[0] =  timesI( Real(metric(mu,tau)) * pp[nu] +  Real(metric(nu,tau)) * pp[mu] )   ;
	out[1] =  timesI( Real(metric(mu,tau)) * pp[nu] -  Real(metric(nu,tau)) * pp[mu] )   ;
	out[2] =  timesI( Real(metric(mu,nu)) * pp[tau] ) ;
	out[3] =  timesI( pp[mu] * pp[nu] * pp[tau] );

	return out;
      };
    int getNumTerms() const
      {
	return 4;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };

  class WvfGamma5NablaMuNablaNuGammaTauA0 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGamma5NablaMuNablaNuGammaTauA0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Real> pp = lower4Vec( make4Vec(mass,p,lattice) );
	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 
	Complex sum = zero;

	for(int alpha=0; alpha < Nd; ++alpha)
	{
	  n1[0] = mu; n1[1] = nu; n1[2] = tau; n1[3] = alpha;
	
	  Real a = peekObs(anti,n1);
	  sum += a * pp[alpha];
	}
	Array<Complex> out(1);
	out[0] = timesI(sum);
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };

  class WvfGamma5NablaMuNablaNuGammaTauPi1 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGamma5NablaMuNablaNuGammaTauPi1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(4);
	Array<Complex> eps = lower4Vec(minkPolVec(mass,p,r,lattice));
	Array<Real> pp = make4Vec(mass,p,lattice);
	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 
	for(int i=0; i < 4; ++i)	{out[i] = zero;}

	for(int alpha=0; alpha < Nd; ++alpha)
	{
	  n1[0] = mu; n1[1] = nu; n1[2] = tau; n1[3] = alpha;
	
	  Real a = peekObs(anti,n1);
	  out[0] += timesI(a * eps[alpha]);
	}
	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int beta=0; beta < Nd; ++beta) 
	  {
	    n1[0] = mu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;
	    Real a = peekObs(anti,n1);
	    Complex d1 = a * eps[beta] * lower4Vec(pp)[alpha] * pp[nu];	    
	   
	    n1[0] = nu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;
	    a = peekObs(anti,n1);
	    Complex d2 = a * eps[beta] * lower4Vec(pp)[alpha] * pp[mu];	 	   
	    
	    out[1] += timesI(d1 + d2);
	    out[2] += timesI(d1 - d2);

	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = beta;
	    a = peekObs(anti,n1);
	    out[3] += timesI(a * eps[beta] * lower4Vec(pp)[alpha] * pp[tau]);
	  }

	return out;
      };
    int getNumTerms() const
      {
	return 4;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfGamma5NablaMuNablaNuGammaTauA1 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGamma5NablaMuNablaNuGammaTauA1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(6);
	Array<Complex> eps = minkPolVec(mass,p,r,lattice);
	Array<Real> pp = make4Vec(mass,p,lattice);

	out[0] = Real(metric(mu,tau)) * eps[nu] +  Real(metric(nu,tau)) * eps[mu];
	out[1] = Real(metric(mu,tau)) * eps[nu] -  Real(metric(nu,tau)) * eps[mu];
	out[2] = Real(metric(mu,nu)) * eps[tau];
	out[3] = pp[mu]* pp[tau] * eps[nu] + pp[nu] * pp[tau] * eps[mu];
	out[4] = pp[mu]* pp[tau] * eps[nu] - pp[nu] * pp[tau] * eps[mu];
	out[5] = pp[mu] * pp[nu] * eps[tau];

	return out;
      };
    int getNumTerms() const
      {
	return 6;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfGamma5NablaMuNablaNuGammaTauPi2 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGamma5NablaMuNablaNuGammaTauPi2(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(3);
	Array< Array<Complex> > tens = minkPolTens(mass,p,r,lattice); 
	Array<Real> pp = make4Vec(mass,p,lattice);

	out[0] = tens[mu][tau] * pp[nu] +  tens[nu][tau] * pp[mu];
	out[1] = tens[mu][tau] * pp[nu] -  tens[nu][tau] * pp[mu];
	out[2] = tens[mu][nu] * pp[tau];

	return out;
      };
    int getNumTerms() const
      {
	return 3;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfGamma5NablaMuNablaNuGammaTauA2 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGamma5NablaMuNablaNuGammaTauA2(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(3);
	for(int i=0; i < 3; ++i)	{out[i] = zero;}

	Array< Array<Complex> > tens = minkPolTens(mass,p,r,lattice); 
	Array<Real> pp = lower4Vec( make4Vec(mass,p,lattice) );

	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 

	Array<Complex> tempmu = lower4Vec( tens[mu] );
	Array<Complex> tempnu = lower4Vec( tens[nu] );
	Array<Complex> temptau = lower4Vec( tens[tau] );

	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int beta=0; beta < Nd; ++beta) 
	  {
	    n1[0] = mu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;	    
	    Real a = peekObs(anti,n1);
	    Complex d1 = a * tempnu[alpha] * pp[beta] ;

	    n1[0] = nu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;	    
	    a = peekObs(anti,n1);
	    Complex d2 = a * tempmu[alpha] * pp[beta] ;
	    
	    out[0] += timesI(d1 + d2);
	    out[1] += timesI(d1 - d2);

	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = beta;	    
	    a = peekObs(anti,n1);
	    out[2] = timesI(a * temptau[alpha] * pp[beta]) ;	    
	  }  

	return out;
      };
    int getNumTerms() const
      {
	return 3;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfGamma5NablaMuNablaNuGammaTauA3 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGamma5NablaMuNablaNuGammaTauA3(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array< Array< Array<Complex> > > tens3 = minkPolTensSpin3(mass,p,r,lattice); 
	Array<Complex> out(1);
	out[0] = tens3[mu][nu][tau];

	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };




  // FMuNu GammaTau


  class WvfFMuNuGammaTauA0 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfFMuNuGammaTauA0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Real> pp = make4Vec(mass,p,lattice);
	Array<Complex> out(4);  
	out[0] =  ( Real(metric(mu,tau)) * pp[nu] +  Real(metric(nu,tau)) * pp[mu] )   ;
	out[1] =  ( Real(metric(mu,tau)) * pp[nu] -  Real(metric(nu,tau)) * pp[mu] )   ;
	out[2] =  Real(metric(mu,nu)) * pp[tau] ;
	out[3] =  pp[mu] * pp[nu] * pp[tau];

	return out;
      };
    int getNumTerms() const
      {
	return 4;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };

  class WvfFMuNuGammaTauPi0 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfFMuNuGammaTauPi0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Real> pp = lower4Vec( make4Vec(mass,p,lattice) );
	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 
	Complex sum = zero;

	for(int alpha=0; alpha < Nd; ++alpha)
	{
	  n1[0] = mu; n1[1] = nu; n1[2] = tau; n1[3] = alpha;
	
	  Real a = peekObs(anti,n1);
	  sum += a * pp[alpha];
	}
	Array<Complex> out(1);
	out[0] = sum;
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };

  class WvfFMuNuGammaTauA1 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfFMuNuGammaTauA1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(4);
	Array<Complex> eps = lower4Vec(minkPolVec(mass,p,r,lattice));
	Array<Real> pp = make4Vec(mass,p,lattice);
	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 
	for(int i=0; i < 4; ++i)	{out[i] = zero;}

	for(int alpha=0; alpha < Nd; ++alpha)
	{
	  n1[0] = mu; n1[1] = nu; n1[2] = tau; n1[3] = alpha;
	
	  Real a = peekObs(anti,n1);
	  out[0] += timesI(a * eps[alpha]);
	}
	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int beta=0; beta < Nd; ++beta) 
	  {
	    n1[0] = mu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;
	    Real a = peekObs(anti,n1);
	    Complex d1 = a * eps[beta] * lower4Vec(pp)[alpha] * pp[nu];	    
	   
	    n1[0] = nu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;
	    a = peekObs(anti,n1);
	    Complex d2 = a * eps[beta] * lower4Vec(pp)[alpha] * pp[mu];	 	   
	    
	    out[1] += timesI(d1 + d2);
	    out[2] += timesI(d1 - d2);

	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = beta;
	    a = peekObs(anti,n1);
	    out[3] += timesI(a * eps[beta] * lower4Vec(pp)[alpha] * pp[tau]);
	  }

	return out;
      };
    int getNumTerms() const
      {
	return 4;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfFMuNuGammaTauPi1 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfFMuNuGammaTauPi1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(6);
	Array<Complex> eps = minkPolVec(mass,p,r,lattice);
	Array<Real> pp = make4Vec(mass,p,lattice);

	out[0] = Real(metric(mu,tau)) * eps[nu] +  Real(metric(nu,tau)) * eps[mu];
	out[1] = Real(metric(mu,tau)) * eps[nu] -  Real(metric(nu,tau)) * eps[mu];
	out[2] = Real(metric(mu,nu)) * eps[tau];
	out[3] = pp[mu]* pp[tau] * eps[nu] + pp[nu] * pp[tau] * eps[mu];
	out[4] = pp[mu]* pp[tau] * eps[nu] - pp[nu] * pp[tau] * eps[mu];
	out[5] = pp[mu] * pp[nu] * eps[tau];

	return out;
      };
    int getNumTerms() const
      {
	return 6;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfFMuNuGammaTauA2 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfFMuNuGammaTauA2(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(3);
	Array< Array<Complex> > tens = minkPolTens(mass,p,r,lattice); 
	Array<Real> pp = make4Vec(mass,p,lattice);

	out[0] = tens[mu][tau] * pp[nu] +  tens[nu][tau] * pp[mu];
	out[1] = tens[mu][tau] * pp[nu] -  tens[nu][tau] * pp[mu];
	out[2] = tens[mu][nu] * pp[tau];

	return out;
      };
    int getNumTerms() const
      {
	return 3;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfFMuNuGammaTauPi2 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfFMuNuGammaTauPi2(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(3);
	for(int i=0; i < 3; ++i)	{out[i] = zero;}

	Array< Array<Complex> > tens = minkPolTens(mass,p,r,lattice); 
	Array<Real> pp = lower4Vec( make4Vec(mass,p,lattice) );

	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 

	Array<Complex> tempmu = lower4Vec( tens[mu] );
	Array<Complex> tempnu = lower4Vec( tens[nu] );
	Array<Complex> temptau = lower4Vec( tens[tau] );

	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int beta=0; beta < Nd; ++beta) 
	  {
	    n1[0] = mu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;	    
	    Real a = peekObs(anti,n1);
	    Complex d1 = a * tempnu[alpha] * pp[beta] ;

	    n1[0] = nu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;	    
	    a = peekObs(anti,n1);
	    Complex d2 = a * tempmu[alpha] * pp[beta] ;
	    
	    out[0] += timesI(d1 + d2);
	    out[1] += timesI(d1 - d2);

	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = beta;	    
	    a = peekObs(anti,n1);
	    out[2] = timesI(a * temptau[alpha] * pp[beta]) ;	    
	  }  

	return out;
      };
    int getNumTerms() const
      {
	return 3;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfFMuNuGammaTauPi3 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfFMuNuGammaTauPi3(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array< Array< Array<Complex> > > tens3 = minkPolTensSpin3(mass,p,r,lattice); 
	Array<Complex> out(1);
	out[0] = tens3[mu][nu][tau];

	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };



  // Gamma5 FMuNu GammaTau


  class WvfGamma5FMuNuGammaTauRho0 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGamma5FMuNuGammaTauRho0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Real> pp = make4Vec(mass,p,lattice);
	Array<Complex> out(4);  
	out[0] =  ( Real(metric(mu,tau)) * pp[nu] +  Real(metric(nu,tau)) * pp[mu] )   ;
	out[1] =  ( Real(metric(mu,tau)) * pp[nu] -  Real(metric(nu,tau)) * pp[mu] )   ;
	out[2] =  Real(metric(mu,nu)) * pp[tau] ;
	out[3] =  pp[mu] * pp[nu] * pp[tau];

	return out;
      };
    int getNumTerms() const
      {
	return 4;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };

  class WvfGamma5FMuNuGammaTauB0 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGamma5FMuNuGammaTauB0(const Real& mass_, const ArrayInt& p_, LatticeParam lattice_) : mass(mass_), p(p_),  lattice(lattice_) {}
    
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Real> pp = lower4Vec( make4Vec(mass,p,lattice) );
	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 
	Complex sum = zero;

	for(int alpha=0; alpha < Nd; ++alpha)
	{
	  n1[0] = mu; n1[1] = nu; n1[2] = tau; n1[3] = alpha;
	
	  Real a = peekObs(anti,n1);
	  sum += a * pp[alpha];
	}
	Array<Complex> out(1);
	out[0] = sum;
	return out;
      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    LatticeParam lattice;
  };

  class WvfGamma5FMuNuGammaTauRho1 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGamma5FMuNuGammaTauRho1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(4);
	Array<Complex> eps = lower4Vec(minkPolVec(mass,p,r,lattice));
	Array<Real> pp =  make4Vec(mass,p,lattice);
	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 
	for(int i=0; i < 4; ++i)	{out[i] = zero;}

	for(int alpha=0; alpha < Nd; ++alpha)
	{
	  n1[0] = mu; n1[1] = nu; n1[2] = tau; n1[3] = alpha;
	
	  Real a = peekObs(anti,n1);
	  out[0] += timesI(a * eps[alpha]);
	}
	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int beta=0; beta < Nd; ++beta) 
	  {
	    n1[0] = mu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;
	    Real a = peekObs(anti,n1);
	    Complex d1 = a * eps[beta] * lower4Vec(pp)[alpha] * pp[nu];	    
	   
	    n1[0] = nu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;
	    a = peekObs(anti,n1);
	    Complex d2 = a * eps[beta] * lower4Vec(pp)[alpha] * pp[mu];	 	   
	    
	    out[1] += timesI(d1 + d2);
	    out[2] += timesI(d1 - d2);

	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = beta;
	    a = peekObs(anti,n1);
	    out[3] += timesI(a * eps[beta] * lower4Vec(pp)[alpha] * pp[tau]);
	  }

	return out;
      };
    int getNumTerms() const
      {
	return 4;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfGamma5FMuNuGammaTauB1 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGamma5FMuNuGammaTauB1(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(6);
	Array<Complex> eps = minkPolVec(mass,p,r,lattice);
	Array<Real> pp = make4Vec(mass,p,lattice);

	out[0] = Real(metric(mu,tau)) * eps[nu] +  Real(metric(nu,tau)) * eps[mu];
	out[1] = Real(metric(mu,tau)) * eps[nu] -  Real(metric(nu,tau)) * eps[mu];
	out[2] = Real(metric(mu,nu)) * eps[tau];
	out[3] = pp[mu]* pp[tau] * eps[nu] + pp[nu] * pp[tau] * eps[mu];
	out[4] = pp[mu]* pp[tau] * eps[nu] - pp[nu] * pp[tau] * eps[mu];
	out[5] = pp[mu] * pp[nu] * eps[tau];

	return out;
      };
    int getNumTerms() const
      {
	return 6;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfGamma5FMuNuGammaTauRho2 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGamma5FMuNuGammaTauRho2(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(3);
	Array< Array<Complex> > tens = minkPolTens(mass,p,r,lattice); 
	Array<Real> pp = make4Vec(mass,p,lattice);

	out[0] = tens[mu][tau] * pp[nu] +  tens[nu][tau] * pp[mu];
	out[1] = tens[mu][tau] * pp[nu] -  tens[nu][tau] * pp[mu];
	out[2] = tens[mu][nu] * pp[tau];

	return out;
      };
    int getNumTerms() const
      {
	return 3;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfGamma5FMuNuGammaTauB2 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGamma5FMuNuGammaTauB2(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array<Complex> out(3);
	for(int i=0; i < 3; ++i)	{out[i] = zero;}

	Array< Array<Complex> > tens = minkPolTens(mass,p,r,lattice); 
	Array<Real> pp = lower4Vec( make4Vec(mass,p,lattice) );

	TensorReal anti = antiSymTensor<4>();
	Array<int> n1(4); 

	Array<Complex> tempmu = lower4Vec( tens[mu] );
	Array<Complex> tempnu = lower4Vec( tens[nu] );
	Array<Complex> temptau = lower4Vec( tens[tau] );

	for(int alpha=0; alpha < Nd; ++alpha)
	  for(int beta=0; beta < Nd; ++beta) 
	  {
	    n1[0] = mu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;	    
	    Real a = peekObs(anti,n1);
	    Complex d1 = a * tempnu[alpha] * pp[beta] ;

	    n1[0] = nu; n1[1] = tau; n1[2] = alpha; n1[3] = beta;	    
	    a = peekObs(anti,n1);
	    Complex d2 = a * tempmu[alpha] * pp[beta] ;
	    
	    out[0] += timesI(d1 + d2);
	    out[1] += timesI(d1 - d2);

	    n1[0] = mu; n1[1] = nu; n1[2] = alpha; n1[3] = beta;	    
	    a = peekObs(anti,n1);
	    out[2] = timesI(a * temptau[alpha] * pp[beta]) ;	    
	  }  

	return out;
      };
    int getNumTerms() const
      {
	return 3;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };

  class WvfGamma5FMuNuGammaTauB3 : public Wvf3index<Complex>
  {
  public:
    //constructor 
    WvfGamma5FMuNuGammaTauB3(const Real& mass_, const ArrayInt& p_, int r_, LatticeParam lattice_) : mass(mass_), p(p_), r(r_),  lattice(lattice_) {}
    
    //symmetrised in mu/nu
    Array<Complex> operator()(int mu, int nu, int tau) const
      {
	Array< Array< Array<Complex> > > tens3 = minkPolTensSpin3(mass,p,r,lattice); 
	Array<Complex> out(1);
	out[0] = tens3[mu][nu][tau];

	return out;

      };
    int getNumTerms() const
      {
	return 1;
      }
  private:
    Real mass;
    ArrayInt p;
    int r;
    LatticeParam lattice;
  };





  //----------------------------------------------------------------------------------
  // Read parameters
  void read(XMLReader& xml, const std::string& path, WaveFuncParams& param)
  {
    WaveFuncParams tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const std::string& path, const WaveFuncParams& param)
  {
    param.writeXML(xml, path);
  }

  //! Initialize
  WaveFuncParams::WaveFuncParams()
  {
  }


  //! Read parameters
  WaveFuncParams::WaveFuncParams(XMLReader& xml, const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "WaveFunctionType",  wavefunc_type);
    read(paramtop, "LatticeParam", lattice);
  }


  // Writer
  void WaveFuncParams::writeXML(XMLWriter& xml, const std::string& path) const
  {
    push(xml, path);

    write(xml, "WaveFunctionType", wavefunc_type);
    write(xml, "LatticeParam", lattice);

    pop(xml);
  }


  //-----------------------------------------------------------------
  //   projection functions for Wvf classes
  //-----------------------------------------------------------------

  Array<Complex> dii(const Wvf2index<Complex>& wvf) 
  {
    int len = wvf.getNumTerms();
    Array<Complex> out(len);
    for(int n=0; n < len; ++n)
    {
      out[n] = zero;
      for(int i=1; i < Nd; ++i) { out[n] += - wvf(i,i)[n]; }   // minus from lowering one spatial index   
    }
    return out;
  }

  Array<Complex> eijk(const Wvf2index<Complex>& wvf, int i)
  {
    int len = wvf.getNumTerms();
    Array<Complex> out(len);
    TensorReal anti = antiSymTensor<3>();
    Array<int> n1(3); 
 
    for(int n=0; n < len; ++n)
    {
      out[n] = zero;
      for(int j=1; j < Nd; ++j)
	for(int k=1; k < Nd; ++k)
	{
	  n1[0] = i - 1 ; n1[1] = j - 1; n1[2] = k - 1;
	      
	  Real a = peekObs(anti,n1);
	  out[n] += a * wvf(j,k)[n];
	}
    }
    return out;  
  }

  Real test_antisym3()
  {
    TensorReal anti = antiSymTensor<3>();
    Array<int> n1(3);
    n1[0] = 0; n1[1] = 1; n1[2] = 2;
    Real a = peekObs(anti,n1);
    return a;
  }

  Array<Complex> sijk(const Wvf2index<Complex>& wvf, int i)
  {
    int len = wvf.getNumTerms();
    Array<Complex> out(len);
    TensorReal anti = antiSymTensor<3>();
    Array<int> n1(3); 
 
    for(int n=0; n < len; ++n)
    {
      out[n] = zero;
      for(int j=1; j < Nd; ++j)
	for(int k=1; k < Nd; ++k)
	{
	  n1[0] = i - 1 ; n1[1] = j - 1; n1[2] = k - 1;
	      
	  Real a = peekObs(anti,n1);
	  out[n] += fabs(a) * wvf(j,k)[n];     //
	}
    }
    return out;  
  }

  Array<Complex> Sajk(const Wvf2index<Complex>& wvf, int a)
  {
    int len = wvf.getNumTerms();
    Array<Complex> out(len);
    TensorReal anti = antiSymTensor<3>();
    Array<int> n1(3); 
 
    for(int n=0; n < len; ++n)
    {
      out[n] = zero;
      for(int j=1; j < Nd; ++j)
	for(int k=1; k < Nd; ++k)
	{
	  Real b = Real(0.0);
	  if (j!=k)
	  {b = Real(0.0);}
	  else if ( ((a==1) && (j==1) && (k==1))  || ((a==2) && (j==2) && (k==2))  )
	  {b = Real(1.0);}
	  else if ( ((a==1) && (j==2) && (k==2))  || ((a==2) && (j==3) && (k==3))  )
	  {b = Real(-1.0);};
	      
	  out[n] += b * wvf(j,k)[n];      
	}
    }
    return out;  
  }

  Array<Complex> d0eijk(const Wvf3index<Complex>& wvf, int dir)
  {
    int len = wvf.getNumTerms();
    Array<Complex> out(len);
    TensorReal anti = antiSymTensor<3>();
    Array<int> n1(3); 

    for(int n=0; n < len; ++n)
    {
      out[n] = zero;

      for(int j=1; j < Nd; ++j)
	for(int k=1; k < Nd; ++k)
	{
	  n1[0] = dir - 1 ; n1[1] = j - 1; n1[2] = k - 1;
	      
	  Real a = peekObs(anti,n1);
	  out[n] += a * wvf(j, k, 0)[n];
	}
    }
    return out;
  }

  Array<Complex> d0sijk(const Wvf3index<Complex>& wvf, int dir)
  {
    int len = wvf.getNumTerms();
    Array<Complex> out(len);
    TensorReal anti = antiSymTensor<3>();
    Array<int> n1(3); 
   
    for(int n=0; n < len; ++n)
    {
      out[n] = zero;
	    
      for(int j=1; j < Nd; ++j)
	for(int k=1; k < Nd; ++k)
	{
	  n1[0] = dir - 1 ; n1[1] = j - 1; n1[2] = k - 1;
	      
	  Real a = peekObs(anti,n1);
	  out[n] += fabs(a) * wvf(j, k, 0)[n];
	}
    }
    return out;
  }
 
  Array<Complex> diieijk(const Wvf3index<Complex>& wvf)
  {
    int len = wvf.getNumTerms();
    Array<Complex> out(len);
    TensorReal anti = antiSymTensor<3>();
    Array<int> n1(3); 

   
    for(int n=0; n < len; ++n)
    {
      out[n] = zero;

      for(int l=1; l < Nd; ++l)
      {
	    
	for(int j=1; j < Nd; ++j)
	  for(int k=1; k < Nd; ++k)
	  {
	    n1[0] = l - 1 ; n1[1] = j - 1; n1[2] = k - 1;
		  
	    Real a = peekObs(anti,n1);
	    out[n] += a * wvf(j, k, l)[n];
	  }
      }
    }
    return out;
  }

  Array<Complex> diisijk(const Wvf3index<Complex>& wvf)
  {
    int len = wvf.getNumTerms();
    Array<Complex> out(len);
    TensorReal anti = antiSymTensor<3>();
    Array<int> n1(3); 

   
    for(int n=0; n < len; ++n)
    {
      out[n] = zero;

      for(int l=1; l < Nd; ++l)
      {
	    
	for(int j=1; j < Nd; ++j)
	  for(int k=1; k < Nd; ++k)
	  {
	    n1[0] = l - 1 ; n1[1] = j - 1; n1[2] = k - 1;
		  
	    Real a = peekObs(anti,n1);
	    out[n] += fabs(a) * wvf(j, k, l)[n];
	  }
      }
    }
    return out;
  }

  Array<Complex> eijkeijk(const Wvf3index<Complex>& wvf, int dir)
  {
    int len = wvf.getNumTerms();
    Array<Complex> out(len);
    TensorReal anti = antiSymTensor<3>();
    Array<int> n1(3); 
    Array<int> n2(3);
   
    for(int n=0; n < len; ++n)
    {
      out[n] = zero;

      for(int l=1; l < Nd; ++l)
	for(int m=1; m < Nd; ++m)
	{
	  n2[0] = dir - 1 ; n2[1] = l - 1 ; n2[2] = m - 1;
	  Real b = peekObs(anti,n2);

	  Complex sum = zero;	      
	  for(int j=1; j < Nd; ++j)
	    for(int k=1; k < Nd; ++k)
	    {
	      n1[0] = m - 1 ; n1[1] = j - 1; n1[2] = k - 1;
		    
	      Real a = peekObs(anti,n1);
	      sum += a * wvf(j, k, l)[n];
	    }
	  out[n] += b * sum;
	}
	
    }
    return out;
  }

  Array<Complex> sijkeijk(const Wvf3index<Complex>& wvf, int dir)
  {
    int len = wvf.getNumTerms();
    Array<Complex> out(len);
    TensorReal anti = antiSymTensor<3>();
    Array<int> n1(3); 
    Array<int> n2(3);
   
    for(int n=0; n < len; ++n)
    {
      out[n] = zero;

      for(int l=1; l < Nd; ++l)
	for(int m=1; m < Nd; ++m)
	{
	  n2[0] = dir - 1 ; n2[1] = l - 1 ; n2[2] = m - 1;
	  Real b = peekObs(anti,n2);

	  Complex sum = zero;	      
	  for(int j=1; j < Nd; ++j)
	    for(int k=1; k < Nd; ++k)
	    {
	      n1[0] = m - 1 ; n1[1] = j - 1; n1[2] = k - 1;
		    
	      Real a = peekObs(anti,n1);
	      sum += a * wvf(j, k, l)[n];
	    }
	  out[n] += fabs(b) * sum;
	}
	
    }
    return out;
  }



  Array<Complex> Sajkeijk(const Wvf3index<Complex>& wvf, int dir)
  {
    int len = wvf.getNumTerms();
    Array<Complex> out(len);
    TensorReal anti = antiSymTensor<3>();
    Array<int> n1(3); 
   
    for(int n=0; n < len; ++n)
    {
      out[n] = zero;

      for(int l=1; l < Nd; ++l)
	for(int m=1; m < Nd; ++m)
	{
	  Real b = Real(0.0);
	  if (l!=m)
	  {b = Real(0.0);}
	  else if ( ((dir==1) && (l==1) && (m==1))  || ((dir==2) && (l==2) && (m==2))  )
	  {b = Real(1.0);}
	  else if ( ((dir==1) && (l==2) && (m==2))  || ((dir==2) && (l==3) && (m==3))  )
	  {b = Real(-1.0);};

	  Complex sum = zero;	      
	  for(int j=1; j < Nd; ++j)
	    for(int k=1; k < Nd; ++k)
	    {
	      n1[0] = m - 1 ; n1[1] = j - 1; n1[2] = k - 1;
		    
	      Real a = peekObs(anti,n1);
	      sum += a * wvf(j, k, l)[n];
	    }
	  out[n] += b * sum;
	}	
    }
    return out;
  }



  Array<Complex> eijksijk(const Wvf3index<Complex>& wvf, int dir)
  {
    int len = wvf.getNumTerms();
    Array<Complex> out(len);
    TensorReal anti = antiSymTensor<3>();
    Array<int> n1(3); 
    Array<int> n2(3);
   
    for(int n=0; n < len; ++n)
    {
      out[n] = zero;

      for(int l=1; l < Nd; ++l)
	for(int m=1; m < Nd; ++m)
	{
	  n2[0] = dir - 1 ; n2[1] = l - 1 ; n2[2] = m - 1;
	  Real b = peekObs(anti,n2);

	  Complex sum = zero;	      
	  for(int j=1; j < Nd; ++j)
	    for(int k=1; k < Nd; ++k)
	    {
	      n1[0] = m - 1 ; n1[1] = j - 1; n1[2] = k - 1;
		    
	      Real a = peekObs(anti,n1);
	      sum += fabs(a) * wvf(j, k, l)[n];
	    }
	  out[n] += b * sum;
	}
	
    }
    return out;
  }

  Array<Complex> sijksijk(const Wvf3index<Complex>& wvf, int dir)
  {
    int len = wvf.getNumTerms();
    Array<Complex> out(len);
    TensorReal anti = antiSymTensor<3>();
    Array<int> n1(3); 
    Array<int> n2(3);
   
    for(int n=0; n < len; ++n)
    {
      out[n] = zero;

      for(int l=1; l < Nd; ++l)
	for(int m=1; m < Nd; ++m)
	{
	  n2[0] = dir - 1 ; n2[1] = l - 1 ; n2[2] = m - 1;
	  Real b = peekObs(anti,n2);
	      
	  Complex sum = zero;
	  for(int j=1; j < Nd; ++j)
	    for(int k=1; k < Nd; ++k)
	    {
	      n1[0] = m - 1 ; n1[1] = j - 1; n1[2] = k - 1;
		    
	      Real a = peekObs(anti,n1);
	      sum += fabs(a) * wvf(j, k, l)[n];
	    }
	  out[n] += fabs(b) * sum;
	}
	
    }
    return out;
  }


  Array<Complex> Sajksijk(const Wvf3index<Complex>& wvf, int dir)
  {
    int len = wvf.getNumTerms();
    Array<Complex> out(len);
    TensorReal anti = antiSymTensor<3>();
    Array<int> n1(3); 
   
    for(int n=0; n < len; ++n)
    {
      out[n] = zero;

      for(int l=1; l < Nd; ++l)
	for(int m=1; m < Nd; ++m)
	{
	  Real b = Real(0.0);
	  if (l!=m)
	  {b = Real(0.0);}
	  else if ( ((dir==1) && (l==1) && (m==1))  || ((dir==2) && (l==2) && (m==2))  )
	  {b = Real(1.0);}
	  else if ( ((dir==1) && (l==2) && (m==2))  || ((dir==2) && (l==3) && (m==3))  )
	  {b = Real(-1.0);};

	  Complex sum = zero;	      
	  for(int j=1; j < Nd; ++j)
	    for(int k=1; k < Nd; ++k)
	    {
	      n1[0] = m - 1 ; n1[1] = j - 1; n1[2] = k - 1;
		    
	      Real a = peekObs(anti,n1);
	      sum += fabs(a) * wvf(j, k, l)[n];
	    }
	  out[n] += b * sum;
	}	
    }
    return out;
  }



  //----------------------------------------------------------------------------------
  // Construct (Gamma(0)) wavefunction overlap onto an  a_0
  Array<Complex> MesGamma0A0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfA0(mass, p, params.lattice)();
  }


  //----------------------------------------------------------------------------------
  // Construct (Gamma(1)) wavefunction overlap onto an  b_0
  Array<Complex> MesGamma1B0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGammaMuB0(mass, p, params.lattice)(1);
  }


  // Construct (Gamma(1)) wavefunction overlap onto an  rho_1
  Array<Complex> MesGamma1Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGammaMuRho1(mass, p, r, params.lattice)(1);
  }

  //----------------------------------------------------------------------------------
  // Construct (Gamma(2)) wavefunction overlap onto an  b_0
  Array<Complex> MesGamma2B0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGammaMuB0(mass, p, params.lattice)(2);
  }


  // Construct (Gamma(2)) wavefunction overlap onto an  rho_1
  Array<Complex> MesGamma2Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGammaMuRho1(mass, p, r, params.lattice)(2);
  }

  //----------------------------------------------------------------------------------
  // Construct (Gamma(3)) wavefunction overlap onto an  b_1
  Array<Complex> MesGamma3B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    Array<Complex> out(1);
    out[0] =  Real(-1.0) * timesI( WvfSigmaMuNuB1(mass, p, r, params.lattice)(1 , 2)[0] ) ; 
    return out;
  }


  // Construct (Gamma(3)) wavefunction overlap onto an  rho_1
  Array<Complex> MesGamma3Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    Array<Complex> out(1);
    out[0] =  Real(-1.0) * timesI( WvfSigmaMuNuRho1(mass, p, r, params.lattice)(1 , 2)[0] ) ; 
    return out;
  }

  //----------------------------------------------------------------------------------
  // Construct (Gamma(4)) wavefunction overlap onto an  b_0
  Array<Complex> MesGamma4B0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGammaMuB0(mass, p, params.lattice)(3);
  }


  // Construct (Gamma(4)) wavefunction overlap onto an  rho_1
  Array<Complex> MesGamma4Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGammaMuRho1(mass, p, r, params.lattice)(3);
  }


  //----------------------------------------------------------------------------------
  // Construct (Gamma(5)) wavefunction overlap onto an  b_1
  Array<Complex> MesGamma5B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    Array<Complex> out(1);
    out[0] =  Real(-1.0) * timesI( WvfSigmaMuNuB1(mass, p, r, params.lattice)(1 , 3)[0] ) ; 
    return out;
  }


  // Construct (Gamma(5)) wavefunction overlap onto an  rho_1
  Array<Complex> MesGamma5Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    Array<Complex> out(1);
    out[0] =  Real(-1.0) * timesI( WvfSigmaMuNuRho1(mass, p, r, params.lattice)(1 , 3)[0] ) ; 
    return out;
  }


  //----------------------------------------------------------------------------------
  // Construct (Gamma(6)) wavefunction overlap onto an  b_1
  Array<Complex> MesGamma6B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    Array<Complex> out(1);
    out[0] =  Real(-1.0) * timesI( WvfSigmaMuNuB1(mass, p, r, params.lattice)(2 , 3)[0] ) ; 
    return out;
  }


  // Construct (Gamma(6)) wavefunction overlap onto an  rho_1
  Array<Complex> MesGamma6Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    Array<Complex> out(1);
    out[0] =  Real(-1.0) * timesI( WvfSigmaMuNuRho1(mass, p, r, params.lattice)(2 , 3)[0] ) ; 
    return out;
  }


  //----------------------------------------------------------------------------------
  // Construct (Gamma(7)) wavefunction overlap onto an  pi_0
  Array<Complex> MesGamma7Pi0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGamma5GammaMuPi0(mass, p, params.lattice)(0);
  }


  // Construct (Gamma(7)) wavefunction overlap onto an  a_1
  Array<Complex> MesGamma7A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGamma5GammaMuA1(mass, p, r, params.lattice)(0);
  }


  //----------------------------------------------------------------------------------
  // Construct (Gamma(8)) wavefunction overlap onto an  b_0
  Array<Complex> MesGamma8B0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGammaMuB0(mass, p, params.lattice)(0);
  }


  // Construct (Gamma(8)) wavefunction overlap onto an  rho_1
  Array<Complex> MesGamma8Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGammaMuRho1(mass, p, r, params.lattice)(0);
  }


  //----------------------------------------------------------------------------------
  // Construct (Gamma(9)) wavefunction overlap onto an  b_1
  Array<Complex> MesGamma9B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    Array<Complex> out(1);
    out[0] =  Real(-1.0) * timesI( WvfSigmaMuNuB1(mass, p, r, params.lattice)(1 , 0)[0] ) ; 
    return out;
  }


  // Construct (Gamma(9)) wavefunction overlap onto an  rho_1
  Array<Complex> MesGamma9Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    Array<Complex> out(1);
    out[0] =  Real(-1.0) * timesI( WvfSigmaMuNuRho1(mass, p, r, params.lattice)(1 , 0)[0] ) ; 
    return out;
  }


  //----------------------------------------------------------------------------------
  // Construct (Gamma(10)) wavefunction overlap onto an  b_1
  Array<Complex> MesGamma10B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    Array<Complex> out(1);
    out[0] =  Real(-1.0) * timesI( WvfSigmaMuNuB1(mass, p, r, params.lattice)(2 , 0)[0] ) ; 
    return out;
  }


  // Construct (Gamma(10)) wavefunction overlap onto an  rho_1
  Array<Complex> MesGamma10Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    Array<Complex> out(1);
    out[0] =  Real(-1.0) * timesI( WvfSigmaMuNuRho1(mass, p, r, params.lattice)(2 , 0)[0] ) ; 
    return out;
  }


  //----------------------------------------------------------------------------------
  // Construct (Gamma(11)) wavefunction overlap onto an  pi_0
  Array<Complex> MesGamma11Pi0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    Array< Complex > out(1);
    out[0] = Real(-1.0) * WvfGamma5GammaMuPi0(mass, p, params.lattice)(3)[0];
    return out;
  }


  // Construct (Gamma(11)) wavefunction overlap onto an  a_1
  Array<Complex> MesGamma11A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    Array< Complex > out(1);
    out[0] = Real(-1.0) * WvfGamma5GammaMuA1(mass, p, r, params.lattice)(3)[0]; 
    return out;
  }


  //----------------------------------------------------------------------------------
  // Construct (Gamma(12)) wavefunction overlap onto an  b_1
  Array<Complex> MesGamma12B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    Array<Complex> out(1);
    out[0] =  Real(-1.0) * timesI( WvfSigmaMuNuB1(mass, p, r, params.lattice)(3 , 0)[0] ) ; 
    return out;
  }


  // Construct (Gamma(12)) wavefunction overlap onto an  rho_1
  Array<Complex> MesGamma12Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    Array<Complex> out(1);
    out[0] =  Real(-1.0) * timesI( WvfSigmaMuNuRho1(mass, p, r, params.lattice)(3 , 0)[0] ) ; 
    return out;
  }


  //----------------------------------------------------------------------------------
  // Construct (Gamma(13)) wavefunction overlap onto an  pi_0
  Array<Complex> MesGamma13Pi0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGamma5GammaMuPi0(mass, p, params.lattice)(2);
  }


  // Construct (Gamma(13)) wavefunction overlap onto an  a_1
  Array<Complex> MesGamma13A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGamma5GammaMuA1(mass, p, r, params.lattice)(2);
  }


  //----------------------------------------------------------------------------------
  // Construct (Gamma(14)) wavefunction overlap onto an  pi_0
  Array<Complex> MesGamma14Pi0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    Array< Complex > out(1);
    out[0] = Real(-1.0) * WvfGamma5GammaMuPi0(mass, p, params.lattice)(1)[0]; 
    return out;
  }


  // Construct (Gamma(14)) wavefunction overlap onto an  a_1
  Array<Complex> MesGamma14A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    Array< Complex > out(1);
    out[0] =  Real(-1.0) * WvfGamma5GammaMuA1(mass, p, r, params.lattice)(1)[0];
    return out;
  }


  //----------------------------------------------------------------------------------
  // Construct (Gamma(15)) wavefunction overlap onto an  pi_0
  Array<Complex> MesGamma15Pi0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGamma5Pi0(mass, p, params.lattice)();
  }




  //looks like these are the names in Manke's table

  //----------------------------------------------------------------------------------
  // Construct (a_0) wavefunction overlap onto an  a_0
  Array<Complex> MesA0A1A0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfA0(mass, p, params.lattice)();
  }


  //----------------------------------------------------------------------------------
  // Construct (pion) wavefunction overlap onto an  pi_0
  Array<Complex> MesPionA1Pi0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGamma5Pi0(mass, p, params.lattice)();
  }


  //----------------------------------------------------------------------------------
  // Construct (rho) wavefunction overlap onto an  b_0
  Array<Complex> MesRhoT1B0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGammaMuB0(mass, p, params.lattice)(dir);
  }


  // Construct (rho) wavefunction overlap onto an  rho_1
  Array<Complex> MesRhoT1Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGammaMuRho1(mass, p, r, params.lattice)(dir);
  }


  //----------------------------------------------------------------------------------
  // Construct (a_1) wavefunction overlap onto an  pi_0
  Array<Complex> MesA1T1Pi0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGamma5GammaMuPi0(mass, p, params.lattice)(dir);
  }

  // Construct (a_1) wavefunction overlap onto an  a_1
  Array<Complex> MesA1T1A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGamma5GammaMuA1(mass, p, r, params.lattice)(dir);
  }


  //----------------------------------------------------------------------------------
  // Construct (b_1) wavefunction overlap onto an  b_1
  Array<Complex> MesB1T1B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    int mu, nu;
    switch (dir)
    {
    case 1:
      mu = 2; nu = 3;    
      break;

    case 2:
      mu = 1; nu = 3;
      break;

    case 3:
      mu = 1; nu = 2;
      break;

    default:
      cerr << __func__ << ": invalid direction" << endl;
      exit(1);
    }

    Array<Complex> out(1);
    out[0] =  Real(-1.0) * timesI( WvfSigmaMuNuB1(mass, p, r, params.lattice)(mu,nu)[0]  );
    return out;
  }

  // Construct (b_1) wavefunction overlap onto an  rho_1
  Array<Complex> MesB1T1Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    int mu, nu;
    switch (dir)
    {
    case 1:
      mu = 2; nu = 3;
      break;

    case 2:
      mu = 1; nu = 3;
      break;

    case 3:
      mu = 1; nu = 2;
      break;

    default:
      cerr << __func__ << ": invalid direction" << endl;
      exit(1);
    }

    Array<Complex> out(1);
    out[0] =  Real(-1.0) * timesI( WvfSigmaMuNuRho1(mass, p, r, params.lattice)(mu,nu)[0]  );
    return out;
  }


  //Manke's derivative based operators

  //----------------------------------------------------------------------------------
  // Construct (PionxNabla_T1) wavefunction overlap onto an  rho_0
  Array<Complex> MesPionxNablaT1Rho0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGamma5NablaMuRho0(mass, p, params.lattice)(dir);
  }


  // Construct (PionxNabla_T1) wavefunction overlap onto a  b_1
  Array<Complex> MesPionxNablaT1B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGamma5NablaMuB1(mass, p, r, params.lattice)(dir);
  }


  //----------------------------------------------------------------------------------
  // Construct (A0xNabla_T1) wavefunction overlap onto a  b_0
  Array<Complex> MesA0xNablaT1B0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfNablaMuB0(mass, p,  params.lattice)(dir);
  }


  // Construct (A0xNabla_T1) wavefunction overlap onto a  rho_1
  Array<Complex> MesA0xNablaT1Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfNablaMuRho1(mass, p, r, params.lattice)(dir);
  }


  //----------------------------------------------------------------------------------
  // Construct (A0_2xNabla_T1) wavefunction overlap onto a  a_0 
  Array<Complex> MesA02xNablaT1A0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGammaMuNablaNuA0(mass, p, params.lattice)(0, dir);
  }


  // Construct (A0_2xNabla_T1) wavefunction overlap onto a  a_1
  Array<Complex> MesA02xNablaT1A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGammaMuNablaNuA1(mass, p, r, params.lattice)(0,dir);
  }


  // Construct (A0_2xNabla_T1) wavefunction overlap onto a  a_2
  Array<Complex> MesA02xNablaT1A2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGammaMuNablaNuA2(mass, p, r, params.lattice)(0, dir);
  }


  // Construct (A0_2xNabla_T1) wavefunction overlap onto a  pi_1
  Array<Complex> MesA02xNablaT1Pi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    return WvfGammaMuNablaNuPi1(mass, p, r, params.lattice)(0,dir);
  }





  //----------------------------------------------------------------------------------
  // Construct (RhoxNabla_A1) wavefunction onto an  a_0
  Array<Complex> MesRhoxNablaA1A0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuNablaNuA0 wvf(mass, p,  params.lattice) ;
    Array<Complex> out = dii(wvf);
    return out;
  }

  // overlap on to an a_1 is trivially zero so we won't code it

  // Construct (RhoxNabla_A1) wavefunction onto an  a_2
  Array<Complex> MesRhoxNablaA1A2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuNablaNuA2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = dii(wvf);
    return out;
  }


  // Construct (RhoxNabla_A1) wavefunction onto an  pi_1
  Array<Complex> MesRhoxNablaA1Pi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuNablaNuPi1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = dii(wvf);
    return out; 
  }


  //--------------------------------------------------------------------------------------------------
  // Construct (RhoxNabla_T1) wavefunction onto an a_1
  Array<Complex> MesRhoxNablaT1A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuNablaNuA1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijk(wvf, dir);
    return out; 
  }


  // Construct (RhoxNabla_T1) wavefunction onto a pi_1
  Array<Complex> MesRhoxNablaT1Pi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuNablaNuPi1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijk(wvf, dir);
    return out; 
  }


  //----------------------------------------------------------------------------------------------
  // Construct (RhoxNabla_T2) wavefunction onto an a_0
  Array<Complex> MesRhoxNablaT2A0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuNablaNuA0 wvf(mass, p, params.lattice) ;
    Array<Complex> out = sijk(wvf, dir);
    return out; 
  }

  // Construct (RhoxNabla_T2) wavefunction onto an a_2
  Array<Complex> MesRhoxNablaT2A2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuNablaNuA2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijk(wvf, dir);
    return out; 
  }

  // Construct (RhoxNabla_T2) wavefunction onto a pi_1
  Array<Complex> MesRhoxNablaT2Pi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuNablaNuPi1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijk(wvf, dir);
    return out; 
  }

  //------------------------------------------------------------------------------------------------
  // Construct (RhoxNabla_E) wavefunction onto an a_0 - NB Manke doesnt use this
  Array<Complex> MesRhoxNablaEA0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuNablaNuA0 wvf(mass, p, params.lattice) ;
    Array<Complex> out = Sajk(wvf, dir);
    return out; 
  }

  // Construct (RhoxNabla_E) wavefunction onto an a_2  - NB Manke doesnt use this
  Array<Complex> MesRhoxNablaEA2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuNablaNuA2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajk(wvf, dir);
    return out; 
  }

  // Construct (RhoxNabla_E) wavefunction onto an pi_1 - NB Manke doesnt use this
  Array<Complex> MesRhoxNablaEPi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuNablaNuPi1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajk(wvf, dir);
    return out; 
  }


  //----------------------------------------------------------------------------------
  // Construct (A1xNabla_A1) wavefunction onto an  rho_0 
  Array<Complex> MesA1xNablaA1Rho0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5GammaMuNablaNuRho0 wvf(mass, p, params.lattice) ;
    Array<Complex> out = dii(wvf);
    return out; 
  }

  // Construct (A1xNabla_A1) wavefunction onto an  rho_2 
  Array<Complex> MesA1xNablaA1Rho2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5GammaMuNablaNuRho2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = dii(wvf);
    return out; 
  }

// Construct (A1xNabla_A1) wavefunction onto an  b_1
  Array<Complex> MesA1xNablaA1B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5GammaMuNablaNuB1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = dii(wvf);
    return out; 
  }

  //-------------------------------------------------------------------------------------------------------
  // Construct (A1xNabla_T1) wavefunction onto an rho_1 - not used by Manke
  Array<Complex> MesA1xNablaT1Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5GammaMuNablaNuRho1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijk(wvf, dir);
    return out; 
  }

  // Construct (A1xNabla_T1) wavefunction onto a b_1 - not used by Manke
  Array<Complex> MesA1xNablaT1B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5GammaMuNablaNuB1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijk(wvf, dir);
    return out; 
  }

  //---------------------------------------------------------------------------
  // Construct (A1xNabla_T2) wavefunction onto an rho_0 
  Array<Complex> MesA1xNablaT2Rho0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5GammaMuNablaNuRho0 wvf(mass, p, params.lattice) ;
    Array<Complex> out = sijk(wvf, dir);
    return out; 
  }

  // Construct (A1xNabla_T2) wavefunction onto an rho_2 
  Array<Complex> MesA1xNablaT2Rho2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5GammaMuNablaNuRho2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijk(wvf, dir);
    return out; 
  }

  // Construct (A1xNabla_T2) wavefunction onto an b1 
  Array<Complex> MesA1xNablaT2B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5GammaMuNablaNuB1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijk(wvf, dir);
    return out; 
  }


  //----------------------------------------------------------------------------------------
  // Construct (A1xNabla_E) wavefunction onto an rho_0
  Array<Complex> MesA1xNablaERho0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5GammaMuNablaNuRho0 wvf(mass, p, params.lattice) ;
    Array<Complex> out = Sajk(wvf, dir);
    return out; 
  }

  // Construct (A1xNabla_E) wavefunction onto an rho_2 
  Array<Complex> MesA1xNablaERho2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5GammaMuNablaNuRho2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajk(wvf, dir);
    return out; 
  }

  // Construct (A1xNabla_E) wavefunction onto an b1
  Array<Complex> MesA1xNablaEB1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5GammaMuNablaNuB1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajk(wvf, dir);
    return out; 
  }


  //--------------------------------------------------------------------------------------------
  // Construct (B1xNabla_A1) wavefunction onto an pi_0
  Array<Complex> MesB1xNablaA1Pi0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuGammaNuNablaTauPi0 wvf(mass, p, params.lattice) ;
    Array<Complex> out = diieijk(wvf);
    return out; 
  }

  // Construct (B1xNabla_A1) wavefunction onto an a_1
  Array<Complex> MesB1xNablaA1A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuGammaNuNablaTauA1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = diieijk(wvf);
    return out; 
  }

  // Construct (B1xNabla_A1) wavefunction onto an pi_2
  Array<Complex> MesB1xNablaA1Pi2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuGammaNuNablaTauPi2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = diieijk(wvf);
    return out; 
  }


  //--------------------------------------------------------------------------------------------
  // Construct (B1xNabla_T1) wavefunction onto an a_0
  Array<Complex> MesB1xNablaT1A0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuGammaNuNablaTauA0 wvf(mass, p, params.lattice) ;
    Array<Complex> out = eijkeijk(wvf, dir);
    return out; 
  }

  // Construct (B1xNabla_T1) wavefunction onto an a_1
  Array<Complex> MesB1xNablaT1A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuGammaNuNablaTauA1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijkeijk(wvf, dir);
    return out; 
  }

  // Construct (B1xNabla_T1) wavefunction onto an pi_1
  Array<Complex> MesB1xNablaT1Pi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuGammaNuNablaTauPi1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijkeijk(wvf, dir);
    return out; 
  }

  // Construct (B1xNabla_T1) wavefunction onto an a_2
  Array<Complex> MesB1xNablaT1A2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuGammaNuNablaTauA2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijkeijk(wvf, dir);
    return out; 
  }

  // Construct (B1xNabla_T1) wavefunction onto an pi_2
  Array<Complex> MesB1xNablaT1Pi2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuGammaNuNablaTauPi2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijkeijk(wvf, dir);
    return out; 
  }


  //--------------------------------------------------------------------------------------------
  // Construct (B1xNabla_T2) wavefunction onto an a_1
  Array<Complex> MesB1xNablaT2A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuGammaNuNablaTauA1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijkeijk(wvf, dir);
    return out; 
  }

  // Construct (B1xNabla_T2) wavefunction onto an pi_1
  Array<Complex> MesB1xNablaT2Pi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuGammaNuNablaTauPi1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijkeijk(wvf, dir);
    return out; 
  }

  // Construct (B1xNabla_T2) wavefunction onto an a_2
  Array<Complex> MesB1xNablaT2A2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuGammaNuNablaTauA2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijkeijk(wvf, dir);
    return out; 
  }

  // Construct (B1xNabla_T2) wavefunction onto an pi_2
  Array<Complex> MesB1xNablaT2Pi2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuGammaNuNablaTauPi2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijkeijk(wvf, dir);
    return out; 
  }

  //--------------------------------------------------------------------------------------------
  // Construct (B1xNabla_E) wavefunction onto an a_1
  Array<Complex> MesB1xNablaEA1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuGammaNuNablaTauA1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajkeijk(wvf, dir);
    return out; 
  }

  // Construct (B1xNabla_E) wavefunction onto an pi_1
  Array<Complex> MesB1xNablaEPi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuGammaNuNablaTauPi1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajkeijk(wvf, dir);
    return out; 
  }

  // Construct (B1xNabla_E) wavefunction onto an a_2
  Array<Complex> MesB1xNablaEA2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuGammaNuNablaTauA2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajkeijk(wvf, dir);
    return out; 
  }

  // Construct (B1xNabla_E) wavefunction onto an pi_2
  Array<Complex> MesB1xNablaEPi2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGammaMuGammaNuNablaTauPi2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajkeijk(wvf, dir);
    return out; 
  }



  //-----------------------------------------------------------------------------------------
  // Construct (A02xD_T2) wavefunction onto an b0
  Array<Complex> MesA02xDT2B0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauB0 wvf(mass, p, params.lattice) ;
    Array<Complex> out = d0sijk(wvf, dir);
    return out; 
  }

  // Construct (A02xD_T2) wavefunction onto a b_1
  Array<Complex> MesA02xDT2B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauB1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = d0sijk(wvf, dir);
    return out; 
  }

  // Construct (A02xD_T2) wavefunction onto a \rho_1 
  Array< Complex> MesA02xDT2Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauRho1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = d0sijk(wvf, dir);
    return out; 
  }

  // Construct (A02xD_T2) wavefunction onto a b_2 
  Array<Complex> MesA02xDT2B2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauB2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = d0sijk(wvf, dir);
    return out; 
  }

  // Construct (A02xD_T2) wavefunction onto a \rho_2 
  Array<Complex> MesA02xDT2Rho2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauRho2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = d0sijk(wvf, dir);
    return out;
  }

  // Construct (A02xD_T2) wavefunction onto a \rho_3 
  Array<Complex> MesA02xDT2Rho3WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauRho3 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = d0sijk(wvf, dir);
    return out;
  }


  //-----------------------------------------------------------------------------------------
  // Construct (rhoxD_A2) wavefunction onto a b_2 
  Array<Complex> MesRhoxDA2B2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauB2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = diisijk(wvf);
    return out;
  }

  // Construct (rhoxD_A2) wavefunction onto a \rho_3 
  Array<Complex> MesRhoxDA2Rho3WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauRho3 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = diisijk(wvf);
    return out;
  }


  //--------------------------------------------------------------------------------
  // Construct (rhoxD_T1) wavefunction onto a b0
  Array<Complex> MesRhoxDT1B0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauB0 wvf(mass, p,  params.lattice) ;
    Array<Complex> out = sijksijk(wvf, dir);
    return out;
  }

  // Construct (rhoxD_T1) wavefunction onto a b1
  Array< Complex> MesRhoxDT1B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauB1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijksijk(wvf, dir);
    return out;
  }

  // Construct (rhoxD_T1) wavefunction onto a \rho1
  Array<Complex> MesRhoxDT1Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauRho1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijksijk(wvf, dir);
    return out;
  }

  // Construct (rhoxD_T1) wavefunction onto a \b2
  Array<Complex> MesRhoxDT1B2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauB2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijksijk(wvf, dir);
    return out;
  }

  // Construct (rhoxD_T1) wavefunction onto a \rho2
  Array<Complex> MesRhoxDT1Rho2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauRho2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijksijk(wvf, dir);
    return out;
  }

  // Construct (rhoxD_T1) wavefunction onto a \rho3
  Array<Complex> MesRhoxDT1Rho3WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauRho3 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijksijk(wvf, dir);
    return out;
  }

  //--------------------------------------------------------------------------------
  // Construct (rhoxD_T2) wavefunction onto a \b0
  Array<Complex> MesRhoxDT2B0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauB0 wvf(mass, p, params.lattice) ;
    Array<Complex> out = eijksijk(wvf, dir);
    return out;
  }

  // Construct (rhoxD_T2) wavefunction onto a \b1
  Array<Complex> MesRhoxDT2B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauB1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijksijk(wvf, dir);
    return out;
  }


  // Construct (rhoxD_T2) wavefunction onto a \rho1
  Array<Complex> MesRhoxDT2Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauRho1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijksijk(wvf, dir);
    return out;
  }

  // Construct (rhoxD_T2) wavefunction onto a \b2
  Array<Complex> MesRhoxDT2B2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauB2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijksijk(wvf, dir);
    return out;
  }

  // Construct (rhoxD_T2) wavefunction onto a \rho2
  Array< Complex> MesRhoxDT2Rho2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauRho2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijksijk(wvf, dir);
    return out;
  }

  // Construct (rhoxD_T2) wavefunction onto a \rho3
  Array< Complex> MesRhoxDT2Rho3WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauRho3 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijksijk(wvf, dir);
    return out;
  }



  //------------------------------------------------------------------------------------------------
  // Construct (rhoxD_E) wavefunction onto a b_1
  Array< Complex> MesRhoxDEB1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauB1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajksijk(wvf, dir);
    return out;
  }

  // Construct (rhoxD_E) wavefunction onto a \rho_1
  Array<Complex> MesRhoxDERho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauRho1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajksijk(wvf, dir);
    return out;
  }

  // Construct (rhoxD_E) wavefunction onto a b_2
  Array<Complex> MesRhoxDEB2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauB2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajksijk(wvf, dir);
    return out;
  }

  // Construct (rhoxD_E) wavefunction onto a \rho_2
  Array<Complex> MesRhoxDERho2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuGammaTauRho2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajksijk(wvf, dir);
    return out;
  }


  //-----------------------------------------------------------------------------------------
  // Construct (a1xD_A2) wavefunction onto a pi_2 
  Array<Complex> MesA1xDA2Pi2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauPi2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = diisijk(wvf);
    return out;
  }

  // Construct (a1xD_A2) wavefunction onto an a_3 
  Array<Complex> MesA1xDA2A3WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauA3 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = diisijk(wvf);
    return out;
  }


  //--------------------------------------------------------------------------------
  // Construct (a1xD_T1) wavefunction onto a \pi_0
  Array<Complex> MesA1xDT1Pi0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauPi0 wvf(mass, p,  params.lattice) ;
    Array<Complex> out = sijksijk(wvf, dir);
    return out;
  }

  // Construct (a1xD_T1) wavefunction onto a \pi1
  Array<Complex> MesA1xDT1Pi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauPi1 wvf(mass, p,r,  params.lattice) ;
    Array<Complex> out = sijksijk(wvf, dir);
    return out;
  }

  // Construct (a1xD_T1) wavefunction onto a a1
  Array<Complex> MesA1xDT1A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauA1 wvf(mass, p,r,  params.lattice) ;
    Array<Complex> out = sijksijk(wvf, dir);
    return out;
  }

  // Construct (a1xD_T1) wavefunction onto a \pi2
  Array<Complex> MesA1xDT1Pi2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauPi2 wvf(mass, p,r,  params.lattice) ;
    Array<Complex> out = sijksijk(wvf, dir);
    return out;
  }

  // Construct (a1xD_T1) wavefunction onto a a2
  Array<Complex> MesA1xDT1A2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauA2 wvf(mass, p, r,  params.lattice) ;
    Array<Complex> out = sijksijk(wvf, dir);
    return out;
  }

  // Construct (a1xD_T1) wavefunction onto a a3
  Array<Complex> MesA1xDT1A3WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauA3 wvf(mass, p, r,  params.lattice) ;
    Array<Complex> out = sijksijk(wvf, dir);
    return out;
  }


  //--------------------------------------------------------------------------------
  // Construct (a1xD_T2) wavefunction onto a \pi0
  Array<Complex> MesA1xDT2Pi0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauPi0 wvf(mass, p,  params.lattice) ;
    Array<Complex> out = eijksijk(wvf, dir);
    return out;
  }

  // Construct (a1xD_T2) wavefunction onto a \pi1
  Array<Complex> MesA1xDT2Pi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauPi1 wvf(mass, p,  r, params.lattice) ;
    Array<Complex> out = eijksijk(wvf, dir);
    return out;
  }

  // Construct (a1xD_T2) wavefunction onto a a1
  Array<Complex> MesA1xDT2A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauA1 wvf(mass, p,  r, params.lattice) ;
    Array<Complex> out = eijksijk(wvf, dir);
    return out;
  }

  // Construct (a1xD_T2) wavefunction onto a \pi2
  Array<Complex> MesA1xDT2Pi2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauPi2 wvf(mass, p,  r, params.lattice) ;
    Array<Complex> out = eijksijk(wvf, dir);
    return out;
  }

  // Construct (a1xD_T2) wavefunction onto a a2
  Array<Complex> MesA1xDT2A2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauA2 wvf(mass, p,  r, params.lattice) ;
    Array<Complex> out = eijksijk(wvf, dir);
    return out;
  }

  // Construct (a1xD_T2) wavefunction onto an a3
  Array<Complex> MesA1xDT2A3WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauA3 wvf(mass, p,  r, params.lattice) ;
    Array<Complex> out = eijksijk(wvf, dir);
    return out;
  }



  //------------------------------------------------------------------------------------------------
  // Construct (a1xD_E) wavefunction onto a \pi_1
  Array<Complex> MesA1xDEPi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauPi1 wvf(mass, p,  r, params.lattice) ;
    Array<Complex> out = Sajksijk(wvf, dir);
    return out;
  }

  // Construct (a1xD_E) wavefunction onto a a_1
  Array<Complex> MesA1xDEA1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauA1 wvf(mass, p,  r, params.lattice) ;
    Array<Complex> out = Sajksijk(wvf, dir);
    return out;
  }

  // Construct (a1xD_E) wavefunction onto a \pi_2
  Array<Complex> MesA1xDEPi2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauPi2 wvf(mass, p,  r, params.lattice) ;
    Array<Complex> out = Sajksijk(wvf, dir);
    return out;
  }

  // Construct (a1xD_E) wavefunction onto a a_2
  Array<Complex> MesA1xDEA2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauA2 wvf(mass, p,  r, params.lattice) ;
    Array<Complex> out = Sajksijk(wvf, dir);
    return out;
  }

  //--------------------------------------------------------------------------------------
  // Construct (pion_2xD_T2) wavefunction onto a pi0
  Array<Complex> MesPion2xDT2Pi0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauPi0 wvf(mass, p,  params.lattice) ;
    Array<Complex> out = d0sijk(wvf, dir);
    return out;
  }

  // Construct (pion_2xD_T2) wavefunction onto a pi1
  Array<Complex> MesPion2xDT2Pi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauPi1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = d0sijk(wvf, dir);
    return out;
  }

  // Construct (pion_2xD_T2) wavefunction onto a a1
  Array<Complex> MesPion2xDT2A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauA1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = d0sijk(wvf, dir);
    return out;
  }

  // Construct (pion_2xD_T2) wavefunction onto a pi2
  Array<Complex> MesPion2xDT2Pi2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauPi2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = d0sijk(wvf, dir);
    return out;
  }

  // Construct (pion_2xD_T2) wavefunction onto a a2
  Array<Complex> MesPion2xDT2A2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauA2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = d0sijk(wvf, dir);
    return out;
  }

  // Construct (pion_2xD_T2) wavefunction onto a a3
  Array<Complex> MesPion2xDT2A3WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuGammaTauA3 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = d0sijk(wvf, dir);
    return out;
  }

  //---------------------------------------------------------------------------------------------------
  // Construct (pionxD_T2) wavefunction onto a pi0
  Array<Complex> MesPionxDT2Pi0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuPi0 wvf(mass, p, params.lattice) ;
    Array<Complex> out = sijk(wvf, dir);
    return out;
  }

  // Construct (pionxD_T2) wavefunction onto a pi2
  Array<Complex> MesPionxDT2Pi2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuPi2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijk(wvf, dir);
    return out;
  }

  // Construct (pionxD_T2) wavefunction onto a a1
  Array<Complex> MesPionxDT2A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5NablaMuNablaNuA1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijk(wvf, dir);
    return out;
  }

  //---------------------------------------------------------------------------------------------------
  // Construct (a0xD_T2) wavefunction onto a a0
  Array<Complex> MesA0xDT2A0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuA0 wvf(mass, p, params.lattice) ;
    Array<Complex> out = sijk(wvf, dir);
    return out;
  }

  // Construct (a0xD_T2) wavefunction onto a a2
  Array<Complex> MesA0xDT2A2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuA2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijk(wvf, dir);
    return out;
  }

  // Construct (a0xD_T2) wavefunction onto a pi1
  Array<Complex> MesA0xDT2Pi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfNablaMuNablaNuPi1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijk(wvf, dir);
    return out;
  }

  //---------------------------------------------------------------------------------------------------
  // Construct (a0xB_T1) wavefunction onto a b1
  Array<Complex> MesA0xBT1B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuB1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijk(wvf, dir);
    return out;
  }

  // Construct (a0xB_T1) wavefunction onto a rho1
  Array<Complex> MesA0xBT1Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuRho1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijk(wvf, dir);
    return out;
  }



  //---------------------------------------------------------------------------------------------------
  // Construct (pionxB_T1) wavefunction onto a rho1
  Array<Complex> MesPionxBT1Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuB1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijk(wvf, dir);
    return out;
  }

  // Construct (pionxB_T1) wavefunction onto a b1
  Array<Complex> MesPionxBT1B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuRho1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijk(wvf, dir);
    return out;
  }


  //---------------------------------------------------------------------------------------------------
  // Construct (a0_2xB_T1) wavefunction onto a pi0
  Array<Complex> MesA02xBT1Pi0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauPi0 wvf(mass, p, params.lattice) ;
    Array<Complex> out = d0eijk(wvf, dir);
    return out;
  }

  // Construct (a0_2xB_T1) wavefunction onto a a1
  Array<Complex> MesA02xBT1A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauA1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = d0eijk(wvf, dir);
    return out;
  }

  // Construct (a0_2xB_T1) wavefunction onto a pi1
  Array<Complex> MesA02xBT1Pi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauPi1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = d0eijk(wvf, dir);
    return out;
  }

  // Construct (a0_2xB_T1) wavefunction onto a a2
  Array<Complex> MesA02xBT1A2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauA2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = d0eijk(wvf, dir);
    return out;
  }

  // Construct (a0_2xB_T1) wavefunction onto a pi2
  Array<Complex> MesA02xBT1Pi2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauPi2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = d0eijk(wvf, dir);
    return out;
  }



  //---------------------------------------------------------------------------------------------------
  // Construct (rhoxB_A1) wavefunction onto a pi0
  Array<Complex> MesRhoxBA1Pi0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauPi0 wvf(mass, p, params.lattice) ;
    Array<Complex> out = diieijk(wvf);
    return out;
  }

  // Construct (rhoxB_A1) wavefunction onto a a1
  Array<Complex> MesRhoxBA1A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauA1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = diieijk(wvf);
    return out;
  }

  // Construct (rhoxB_A1) wavefunction onto a pi2
  Array<Complex> MesRhoxBA1Pi2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauPi2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = diieijk(wvf);
    return out;
  }


  //---------------------------------------------------------------------------------------------------
  // Construct (rhoxB_T1) wavefunction onto a a0
  Array<Complex> MesRhoxBT1A0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauA0 wvf(mass, p, params.lattice) ;
    Array<Complex> out = eijkeijk(wvf, dir);
    return out;
  }

  // Construct (rhoxB_T1) wavefunction onto a a1
  Array<Complex> MesRhoxBT1A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauA1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijkeijk(wvf, dir);
    return out;
  }

  // Construct (rhoxB_T1) wavefunction onto a pi1
  Array<Complex> MesRhoxBT1Pi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauPi1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijkeijk(wvf, dir);
    return out;
  }

  // Construct (rhoxB_T1) wavefunction onto a a2
  Array<Complex> MesRhoxBT1A2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauA2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijkeijk(wvf, dir);
    return out;
  }

  // Construct (rhoxB_T1) wavefunction onto a pi2
  Array<Complex> MesRhoxBT1Pi2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauPi2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijkeijk(wvf, dir);
    return out;
  }


  //---------------------------------------------------------------------------------------------------
  // Construct (rhoxB_T2) wavefunction onto a a1
  Array<Complex> MesRhoxBT2A1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauA1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijkeijk(wvf, dir);
    return out;
  }

  // Construct (rhoxB_T2) wavefunction onto a pi1
  Array<Complex> MesRhoxBT2Pi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauPi1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijkeijk(wvf, dir);
    return out;
  }

  // Construct (rhoxB_T2) wavefunction onto a a2
  Array<Complex> MesRhoxBT2A2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauA2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijkeijk(wvf, dir);
    return out;
  }

  // Construct (rhoxB_T2) wavefunction onto a pi2
  Array<Complex> MesRhoxBT2Pi2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauPi2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijkeijk(wvf, dir);
    return out;
  }


  //---------------------------------------------------------------------------------------------------
  // Construct (rhoxB_E) wavefunction onto a a1
  Array<Complex> MesRhoxBEA1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauA1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajkeijk(wvf, dir);
    return out;
  }

  // Construct (rhoxB_E) wavefunction onto a pi1
  Array<Complex> MesRhoxBEPi1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauPi1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajkeijk(wvf, dir);
    return out;
  }

  // Construct (rhoxB_E) wavefunction onto a a2
  Array<Complex> MesRhoxBEA2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauA2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajkeijk(wvf, dir);
    return out;
  }

  // Construct (rhoxB_E) wavefunction onto a pi2
  Array<Complex> MesRhoxBEPi2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfFMuNuGammaTauPi2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajkeijk(wvf, dir);
    return out;
  }


 
  //---------------------------------------------------------------------------------------------------
  // Construct (a1xB_A1) wavefunction onto a b0
  Array<Complex> MesA1xBA1B0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuGammaTauB0 wvf(mass, p, params.lattice) ;
    Array<Complex> out = diieijk(wvf);
    return out;
  }

  // Construct (a1xB_A1) wavefunction onto a rho1
  Array<Complex> MesA1xBA1Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuGammaTauRho1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = diieijk(wvf);
    return out;
  }

  // Construct (a1xB_A1) wavefunction onto a b2
  Array<Complex> MesA1xBA1B2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuGammaTauB2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = diieijk(wvf);
    return out;
  }


  //---------------------------------------------------------------------------------------------------
  // Construct (a1xB_T1) wavefunction onto a rho0
  Array<Complex> MesA1xBT1Rho0WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuGammaTauRho0 wvf(mass, p, params.lattice) ;
    Array<Complex> out = eijkeijk(wvf, dir);
    return out;
  }

  // Construct (a1xB_T1) wavefunction onto a rho1
  Array<Complex> MesA1xBT1Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuGammaTauRho1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijkeijk(wvf, dir);
    return out;
  }

  // Construct (a1xB_T1) wavefunction onto a b1
  Array<Complex> MesA1xBT1B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuGammaTauB1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijkeijk(wvf, dir);
    return out;
  }

  // Construct (a1xB_T1) wavefunction onto a rho2
  Array<Complex> MesA1xBT1Rho2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuGammaTauRho2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijkeijk(wvf, dir);
    return out;
  }

  // Construct (a1xB_T1) wavefunction onto a b2
  Array<Complex> MesA1xBT1B2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuGammaTauB2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = eijkeijk(wvf, dir);
    return out;
  }


  //---------------------------------------------------------------------------------------------------
  // Construct (a1xB_T2) wavefunction onto a rho1
  Array<Complex> MesA1xBT2Rho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuGammaTauRho1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijkeijk(wvf, dir);
    return out;
  }

  // Construct (a1xB_T2) wavefunction onto a b1
  Array<Complex> MesA1xBT2B1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuGammaTauB1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijkeijk(wvf, dir);
    return out;
  }

  // Construct (a1xB_T2) wavefunction onto a rho2
  Array<Complex> MesA1xBT2Rho2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuGammaTauRho2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijkeijk(wvf, dir);
    return out;
  }

  // Construct (a1xB_T2) wavefunction onto a b2
  Array<Complex> MesA1xBT2B2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuGammaTauB2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = sijkeijk(wvf, dir);
    return out;
  }


  //---------------------------------------------------------------------------------------------------
  // Construct (a1xB_E) wavefunction onto a rho1
  Array<Complex> MesA1xBERho1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuGammaTauRho1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajkeijk(wvf, dir);
    return out;
  }

  // Construct (a1xB_E) wavefunction onto a b1
  Array<Complex> MesA1xBEB1WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuGammaTauB1 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajkeijk(wvf, dir);
    return out;
  }

  // Construct (a1xB_E) wavefunction onto a rho2
  Array<Complex> MesA1xBERho2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuGammaTauRho2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajkeijk(wvf, dir);
    return out;
  }

  // Construct (a1xB_E) wavefunction onto a b2
  Array<Complex> MesA1xBEB2WaveFunc::operator()(const Real& mass, const ArrayInt& p, int dir, int r) const
  {
    WvfGamma5FMuNuGammaTauB2 wvf(mass, p, r, params.lattice) ;
    Array<Complex> out = Sajkeijk(wvf, dir);
    return out;
  }



} // namespace FF
