// $Id: formfac_qsq.cc,v 2.0 2008/12/05 04:43:36 edwards Exp $

#include <map>
#include "formfac/formfac_qsq.h"

namespace FF
{

  //! Anonymous namespace
  namespace
  {
    std::ostream& operator<<(std::ostream& s, const ArrayInt& d)
    {
      if (d.size() > 0)
      {
	s << d[0];
	for(int i=1; i < d.size(); ++i)
	  s << " " << d[i];
      }

      return s;
    }


    std::ostream& operator<<(std::ostream& s, const ArrayDouble& d)
    {
      if (d.size() > 0)
      {
	s << d[0];
	for(int i=1; i < d.size(); ++i)
	  s << " " << d[i];
      }

      return s;
    }
  }



  //! 4D dot product
  EnsemReal dot(const EnsemReal& E1, const ArrayDouble& p1, 
		const EnsemReal& E2, const ArrayDouble& p2)
  {
    double sum = 0.0;
    for(int i=0; i < p1.size(); ++i)
      sum += p1[i] * p2[i];

    return (E1*E2 - Real(sum));
  }


  //! 4D Minkowski q^2
  EnsemReal minkQsq(const EnsemReal& E1, const ArrayDouble& p1, 
		    const EnsemReal& E2, const ArrayDouble& p2)
  {
    return pow(E1 - E2,2) - Real(norm2(p1 - p2));
  }


  //! 4D Minkowski q^2 using continuum disp. relation in physical units
  double qsqContDispPhys(double m1, const ArrayInt& p1, 
			 double m2, const ArrayInt& p2,
			 const LatticeParam& lattice)
  {
    const int ns           =  lattice.latt_size[0];   // assumes decay_dir != 0
    const double hc        =  0.197327053; // in GeV fm (conversion factor)
    const double a_s_inv   =  hc / lattice.a_fm;
    const double asq_inv   =  a_s_inv * a_s_inv;
    const double c_sq      =  lattice.aniso.c_sq;
    const double xi        =  lattice.aniso.xi;
    const double xi_sq     =  xi * xi;

    double E1_c  = sqrt( (c_sq/xi_sq)*norm2(contMom(p1,ns)) + m1*m1  );
    double E2_c  = sqrt( (c_sq/xi_sq)*norm2(contMom(p2,ns)) + m2*m2  );

    double q = asq_inv * (xi_sq*(E1_c - E2_c)*(E1_c - E2_c) - norm2(contMom(p1-p2,ns)));

    return q;
  }

  //! 4D Minkowski q^2 using lattice disp. relation in physical units
  double qsqLattDispPhys(double m1, const ArrayInt& p1, 
			 double m2, const ArrayInt& p2,
			 const LatticeParam& lattice)
  {
    const int ns           =  lattice.latt_size[0];   // assumes decay_dir != 0
    const double hc        =  0.197327053; // in GeV fm (conversion factor)
    const double a_s_inv   =  hc / lattice.a_fm;
    const double asq_inv   =  a_s_inv * a_s_inv;
    const double c_sq      =  lattice.aniso.c_sq;
    const double xi        =  lattice.aniso.xi;
    const double xi_sq     =  xi * xi;

    // NOTE: the lattice disp is not really correct for aniso.
    // Need to fix via  anisoHQ.tex  notes
    double E1_l  = acosh( (c_sq/xi_sq)*0.5*norm2(lattMom(p1,ns)) + cosh( m1 ) );
    double E2_l  = acosh( (c_sq/xi_sq)*0.5*norm2(lattMom(p2,ns)) + cosh( m2 ) );

    return asq_inv * (xi_sq*(E1_l - E2_l)*(E1_l - E2_l) - norm2(lattMom(p1-p2,ns)));
  }


  //! Make a continuum momentum from integers
  ArrayDouble contMom(const ArrayInt& nf, int ns)
  {
    ArrayDouble p_c(nf.size());
    const double pi =  3.14159265359;
    
    for(int i=0; i < nf.size(); ++i)
      p_c[i] = 2.0 * pi * nf[i] / ns ;

    return p_c;
  }
  

  //! Make a lattice momentum from integers
  ArrayDouble lattMom(const ArrayInt& nf, int ns)
  {
    ArrayDouble p_l(nf.size());
    const double pi =  3.14159265359;
    
    for(int i=0; i < nf.size(); ++i)
      p_l[i] = 2.0 * sin( pi * nf[i] / ns );

    return p_l;
  }


  //! Canonically order an array of momenta
  /*! \return abs(mom[0]) >= abs(mom[1]) >= ... >= abs(mom[mu]) >= ... >= 0 */
  ArrayInt canonicalOrder(const ArrayInt& mom)
  {
    // first step: make all the components positive
    ArrayInt mom_tmp = mom;
    for (int mu=0; mu < mom_tmp.size(); ++mu)
      if (mom_tmp[mu] < 0) mom_tmp[mu] = -mom_tmp[mu];

    // Initially, the first item is considered sorted.  mu divides mom
    // into a sorted region (<mu) and an unsorted one (>=mu)
    for (int mu=1; mu < mom_tmp.size(); ++mu) 
    {
      // Select the item at the beginning of the unsorted region
      int v = mom_tmp[mu];
      // Work backwards, finding where v should go
      int nu = mu;
      // If this element is less than v, move it up one
      while (mom_tmp[nu-1] < v) {
	mom_tmp[nu] = mom_tmp[nu-1];
	--nu;
	if (nu < 1) break;
      }
      // Stopped when mom_tmp[nu-1] >= v, so put v at postion nu
      mom_tmp[nu] = v;
    }

    return mom_tmp;
  }


  //! List pretty printer
  std::ostream& operator<<(std::ostream& s, const std::list<QsqVal>& qsq_list)
  {
    int j = 0;
    for(std::list<QsqVal>::const_iterator mm = qsq_list.begin();
	mm != qsq_list.end();
	++mm)
    {
      const std::list<PiPf>& pi_pf =  mm->pi_pf;

      s << "int= " << j++
	<< "\tQsq= " << mm->Qsq_c
	<< "\tnum_mom= " << mm->pi_pf.size() 
	<< "\tq_sq= "  << norm2(pi_pf.front().p_f  -  pi_pf.front().p_i)
	<< "\n";

//	<< "    pf_sq= " << norm2(pi_pf.front().p_f)
//	<< "    pi_sq= " << norm2(pi_pf.front().p_i)

      for(std::list<PiPf>::const_iterator nqq=pi_pf.begin(); 
	  nqq != pi_pf.end(); 
	  ++nqq)
      {
	s << "\tp_f= " << nqq->p_f 
	  << "\tp_i= " << nqq->p_i 
	  << std::endl;
      }
    }

    return s;
  }


  //! Writer
  void write(XMLWriter& xml, const std::string& path, const QsqVal& mm)
  {
    push(xml, path);

    write(xml, "Qsq", mm.Qsq_c);
    write(xml, "Qsq_c", mm.Qsq_c);
    write(xml, "Qsq_l", mm.Qsq_l);
    write(xml, "q_sq", norm2(mm.pi_pf.front().p_f  -  mm.pi_pf.front().p_i));

    push(xml, "PiPf");
    for(std::list<PiPf>::const_iterator nqq=mm.pi_pf.begin(); 
	nqq != mm.pi_pf.end(); 
	++nqq)
    {
      write(xml, "elem", *nqq);
    }

    pop(xml);  // PiPf

    pop(xml);  // path
  }


  //! Construct the map of all Qsq points
  std::list<QsqVal> constructQsq(const QsqParam& param)
  {
    const long int int_factor = 100000000;
    int nq_norm_sq_max  = param.nq_norm_sq_max;
    int ni_large = param.ni_norm_sq_max;
    const double mass_i = param.mass_i;
    const double mass_f = param.mass_f;

    //##############################################################################
    //# Figure out ni_max
    //##############################################################################
    int ni_norm_sq_max = nq_norm_sq_max;

    for(std::list<ArrayInt>::const_iterator nf=param.nf_list.begin(); 
	nf != param.nf_list.end(); 
	++nf)
    {
      int nf_norm_sq = norm2(*nf);

      if (ni_norm_sq_max < nq_norm_sq_max + nf_norm_sq)
	ni_norm_sq_max = nq_norm_sq_max + nf_norm_sq;
    }

    //cout << "ni_norm_sq_max = " << ni_norm_sq_max << std::endl;
    int ni_max = int(sqrt(double(ni_norm_sq_max))) + 1;

    int mom_vol = 1;
    ArrayInt mom_size(Nd-1);
    ArrayInt ni_origin(Nd-1);
    ni_origin = -ni_max;
    for(int mu=0; mu < mom_size.size(); ++mu) 
    {
      mom_vol      *= 2*ni_max + 1;
      mom_size[mu]  = 2*ni_max + 1;
    }

    //############################################################################
    // This is totally inefficient but I don't care, it runs fast enough
    //############################################################################

    std::map<long int,QsqVal>  qsq_map;
  
    for(std::list<ArrayInt>::const_iterator nff=param.nf_list.begin(); 
	nff != param.nf_list.end(); 
	++nff)
    {
      const ArrayInt& nf = *nff;

// cout << "nf= " << nf << std::endl;

      for(int n=0; n < mom_vol; ++n)
      {
	ArrayInt ni = crtesn(n, mom_size) + ni_origin;

	// is norm(pf - pi) > norm(q_max) ?

	ArrayInt nq = nf - ni;
	int nq_norm_sq = norm2(nq);

	if ( nq_norm_sq > nq_norm_sq_max ) 
	  continue;

	//check the pi isn't too large
	if( norm2(ni) > ni_large)
	  continue;
	

// cout << "nf= " << nf << "   ni= " << ni << "   nq= " << nq << std::endl;

	// compute Qsq_c
	double Qsq_c = -qsqContDispPhys(mass_f,nf, mass_i,ni, param.lattice);
	double Qsq_l = -qsqLattDispPhys(mass_f,nf, mass_i,ni, param.lattice);

	long int Qsq_c_int = long(int_factor * Qsq_c);
	long int Qsq_l_int = long(int_factor * Qsq_l);

	long int key = Qsq_c_int + Qsq_l_int;  // hopefully this is a unique enough key
//	cout << "Qsq_l= " << Qsq_l << "   Qsq_c= " << Qsq_c << std::endl;
//	cout << "Qsq_l_int= " << Qsq_l_int << "   Qsq_c_int= " << Qsq_c_int << std::endl;
//	cout << "key= " << key << std::endl;

	std::map<long int,QsqVal>::iterator mm = qsq_map.find(key);
	if (mm == qsq_map.end())
	{
	  // This is a new entry
	  QsqVal        q2;
	  PiPf  pi_pf;
	  pi_pf.p_i = ni;
	  pi_pf.p_f = nf;

	  q2.pi_pf.push_back(pi_pf);
	  q2.qsq    = nq_norm_sq;
	  q2.Qsq_c  = Qsq_c;
	  q2.Qsq_l  = Qsq_l;

	  qsq_map[key] = q2;  // insert
	}
	else
	{
	  // Add to a current entry
	  PiPf  pi_pf;
	  pi_pf.p_i = ni;
	  pi_pf.p_f = nf;

	  // Sanity check
	  if (nq_norm_sq != mm->second.qsq)
	  {
	    std::cerr << __func__ 
		      << ": internal error - found an inconsistent 3-vector norm" 
		      << std::endl;
	    exit(1);
	  }
	  
	  mm->second.pi_pf.push_back(pi_pf);
	}
      } // for n
    } // for list nff


    // Return a list
    // NOTE: I should probably sort on something like Qsq_c to be compatible with GTF
    std::list<QsqVal> qsq_list;

    for(std::map<long int,QsqVal>::const_iterator mm = qsq_map.begin();
	mm != qsq_map.end();
	++mm)
    {
      qsq_list.push_back(mm->second);
    }

    return qsq_list;
  }




  //! Momentum norms of source, sink and insertion
  struct FixedMom
  {
    ArrayInt  q;
    ArrayInt  p_i;
    ArrayInt  p_f;
  };


  //! Support for maps
  bool operator<(const FixedMom& a, const FixedMom& b)
  {
    return concat(concat(a.q, a.p_i), a.p_f) < concat(concat(b.q, b.p_i), b.p_f);
  }


  //! Construct the map of all Qsq points
  /* 
   * Here sink, source and insertion 3 momenta are all constrained to fixed
   * norms
   */
  std::list<QsqVal> constructQsqFixedMom(const QsqParam& param)
  {
    const double mass_i    =  param.mass_i;
    const double mass_f    =  param.mass_f;
    const int    nq_norm_sq_max = param.nq_norm_sq_max;

    //##############################################################################
    //# Figure out ni_max
    //##############################################################################
    int ni_norm_sq_max = nq_norm_sq_max;

    for(std::list<ArrayInt>::const_iterator nf=param.nf_list.begin(); 
	nf != param.nf_list.end(); 
	++nf)
    {
      int nf_norm_sq = norm2(*nf);

      if (ni_norm_sq_max < nq_norm_sq_max + nf_norm_sq)
	ni_norm_sq_max = nq_norm_sq_max + nf_norm_sq;
    }

    //cout << "ni_norm_sq_max = " << ni_norm_sq_max << std::endl;
    int ni_max = int(sqrt(double(ni_norm_sq_max))) + 1;

    int mom_vol = 1;
    ArrayInt mom_size(Nd-1);
    ArrayInt ni_origin(Nd-1);
    ni_origin = -ni_max;
    for(int mu=0; mu < mom_size.size(); ++mu) 
    {
      mom_vol      *= 2*ni_max + 1;
      mom_size[mu]  = 2*ni_max + 1;
    }

    //############################################################################
    // This is totally inefficient but I don't care, it runs fast enough
    //############################################################################

    std::map<FixedMom,QsqVal>  qsq_map;
  
    for(std::list<ArrayInt>::const_iterator nff=param.nf_list.begin(); 
	nff != param.nf_list.end(); 
	++nff)
    {
      const ArrayInt& nf = *nff;

// cout << "nf= " << nf << std::endl;

      for(int n=0; n < mom_vol; ++n)
      {
	ArrayInt ni = crtesn(n, mom_size) + ni_origin;

	// is norm(pf - pi) > norm(q_max) ?

	ArrayInt nq = nf - ni;
	int nq_norm_sq = norm2(nq);

	if ( nq_norm_sq > nq_norm_sq_max ) 
	  continue;

// cout << "nf= " << nf << "   ni= " << ni << "   nq= " << nq << std::endl;

	FixedMom key;
	key.q   = canonicalOrder(nq);
	key.p_i = canonicalOrder(ni);
	key.p_f = canonicalOrder(nf);

	std::map<FixedMom,QsqVal>::iterator mm = qsq_map.find(key);
	if (mm == qsq_map.end())
	{
	  // compute Qsq_c
	  double Qsq_c = -qsqContDispPhys(mass_f,nf, mass_i,ni, param.lattice);
	  double Qsq_l = -qsqLattDispPhys(mass_f,nf, mass_i,ni, param.lattice);

//	  std::cout << "nf= " << nf << "  ni= " << ni << "  nq= " << nq << "  Ef=" << Ef_c 
//		    << "  Ei=" << Ei_c << "  Qsq_c =" << Qsq_c << std::endl;

	  // This is a new entry
	  QsqVal  q2;
	  PiPf    pi_pf;
	  pi_pf.p_i = ni;
	  pi_pf.p_f = nf;

	  q2.pi_pf.push_back(pi_pf);
	  q2.qsq    = nq_norm_sq;
	  q2.Qsq_c  = Qsq_c;
	  q2.Qsq_l  = Qsq_l;

	  qsq_map[key] = q2;  // insert
	}
	else
	{
	  // Add to a current entry
	  PiPf  pi_pf;
	  pi_pf.p_i = ni;
	  pi_pf.p_f = nf;

	  // Sanity check
	  if (nq_norm_sq != mm->second.qsq)
	  {
	    std::cerr << __func__ 
		      << ": internal error - found an inconsistent 3-vector norm" 
		      << std::endl;
	    exit(1);
	  }
	  
	  mm->second.pi_pf.push_back(pi_pf);
	}
      } // for n
    } // for list nff


    // Return a list
    // NOTE: I should probably sort on something like Qsq_c to be compatible with GTF
    std::list<QsqVal> qsq_list;

    for(std::map<FixedMom,QsqVal>::const_iterator mm = qsq_map.begin();
	mm != qsq_map.end();
	++mm)
    {
      qsq_list.push_back(mm->second);
    }

    return qsq_list;
  }


} // namespace FF
