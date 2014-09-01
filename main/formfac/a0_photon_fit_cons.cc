 // $Id: a0_photon_fit_cons.cc,v 2.0 2008/12/05 04:43:47 edwards Exp $
// Vector-pseudoscalar form-factor code using fitting method

#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "formfac/formfac.h"
#include "state.h"
#include <list>
#include <iostream>

#include <sys/stat.h>
#include <sys/types.h>


namespace FF
{
  //! Make a state
  class A0Photon3PtFunc : public State3PtFunc
  {
  public:
    //! Constructor
    A0Photon3PtFunc(const std::string& pattern_) : pattern(pattern_) {}

    //! Construct 3pt names
    std::string operator()(int cfg, const ThreePtArg& arg) const
      {
	ArrayInt p_i = arg.pi_pf.p_i;
	ArrayInt p_f = arg.pi_pf.p_f;
	ArrayInt q   = p_f - p_i;
  
	// cout << __func__ << ": pattern=" << pattern.c_str() << std::endl;

	char qxsign = (q[0] < 0 ) ? '-' : '+';
	char qysign = (q[1] < 0 ) ? '-' : '+';
	char qzsign = (q[2] < 0 ) ? '-' : '+';

	char pfxsign = (p_f[0] < 0 ) ? '-' : '+';
	char pfysign = (p_f[1] < 0 ) ? '-' : '+';
	char pfzsign = (p_f[2] < 0 ) ? '-' : '+';

	if (arg.snk[0] < 0 || arg.snk[0] >= 4)
	{
	  std::cerr << __func__ << ": src out of bounds\n";
	  exit(1);
	}
	char xyz[4] = {'X', 'Y', 'Z', 'T'};

	char filen[1024];
	sprintf(filen, pattern.c_str(),
		pfzsign, abs(p_f[2]), pfysign, abs(p_f[1]), pfxsign, abs(p_f[0]),
		cfg,
		xyz[arg.snk[0]],
		qzsign, abs(q[2]), qysign, abs(q[1]), qxsign, abs(q[0]));

	std::cout << __func__ << ": filen=" << filen << "\n";

	return std::string(filen);
      }

  private:
    std::string pattern;
  };

} // namespace FF



//using namespace ADATXML;
using namespace ENSEM;
using namespace FF;
using namespace std;



//! Print an array
std::ostream& operator<<(std::ostream& s, const ArrayInt& d)
{
  s << d[0];
  for(int i=1; i < d.size(); ++i)
    s << " " << d[i];

  return s;
}


//! Print an array
std::ostream& operator<<(std::ostream& s, const ArrayDouble& d)
{
  s << d[0];
  for(int i=1; i < d.size(); ++i)
    s << " " << d[i];

  return s;
}


//! Params for linear system
struct A0PhotonFFSysParam
{
  A0PhotonFFSysParam();
  A0PhotonFFSysParam(XMLReader& xml_in, const std::string& path);

  struct Manage3PtParam
  {
    std::string  pattern;
    std::string  cache_file;
    std::string  cfg_file;
  } threept;

  struct Manage2PtParam
  {
    //  Array<std::string>   pattern_Z_photon;
    Array<std::string>   pattern_Z_a0;
    bool                 avg_2pt_func;
  } twopt;

  struct ManageEnergyParam
  {
    //    Array<std::string>   pattern_photon;
    Real                 photon_Qsq;
    Array<std::string>   pattern_a0;
  } E;

  LatticeParam       lattice;

  string    Z_V_file_name;
  string    qsq_file_name;

  int       t_source;
  int       t_sink;

  struct MomentumRange
  {
    MomentumRange() {mom_rangeP = false; mom2_max = -1;}

    bool             mom_rangeP;
    int              mom2_max;
    Array<ArrayInt>  sink_mom;
  } mom_range;

  struct FixedMomentum
  {
    FixedMomentum() {fixed_momP = false;}

    bool             fixed_momP;
    Array<PiPf>      pi_pf;
  } fixed_mom;

  Array<ArrayInt>  fit_range;

  struct SaveSolutions
  {
    SaveSolutions() {saveP = false;}

    bool          saveP;
    std::string   dir;
    std::string   pattern;
  } saveSoln;


};



#if 0
// NOT READY YET

//! Structure holding mapping of index "unknown" to physical quantities
struct IndexToPhysical
{
  IndexToPhysical() {source_level = sink_level = 0;}

  string       name;  // name of quantity
  int          source_level;
  int          sink_level;

  enum Direction_t
  {
    SPACE_DIRECTION;
    TIME_DIRECTION;
  };

  Direction_t  direction;
};
#endif





// Photon-pi sys params
A0PhotonFFSysParam::A0PhotonFFSysParam(XMLReader& xml_in, const std::string& path) 
{
  try 
  {
    XMLReader paramtop(xml_in, path);

    {
      XMLReader top(paramtop, "ThreePt");

      read(top, "pattern", threept.pattern);
      read(top, "cache_file", threept.cache_file);
      read(top, "cfg_file", threept.cfg_file);
    }

    {
      XMLReader top(paramtop, "TwoPt");

      // read(top, "pattern_Z_photon", twopt.pattern_Z_rho);
      read(top, "pattern_Z_a0", twopt.pattern_Z_a0);
      read(top, "avg_2pt_func", twopt.avg_2pt_func);
    }

    {
      XMLReader top(paramtop, "Energy");

      read(top, "LSZ_photon_Qsq", E.photon_Qsq);
      read(top, "pattern_a0", E.pattern_a0);
    }

    read(paramtop, "LatticeParam", lattice);
    read(paramtop, "Z_V_file_name", Z_V_file_name);
    read(paramtop, "qsq_file_name", qsq_file_name);
    read(paramtop, "t_source", t_source);
    read(paramtop, "t_sink", t_sink);
    read(paramtop, "fit_range", fit_range);

    if (paramtop.count("MomentumRange") == 1)
    {
      XMLReader top(paramtop, "MomentumRange");

      mom_range.mom_rangeP = true;
      read(top, "mom2_max", mom_range.mom2_max);
      read(top, "sink_mom", mom_range.sink_mom);
    }

    if (paramtop.count("FixedMomentum") == 1)
    {
      XMLReader top(paramtop, "FixedMomentum");

      fixed_mom.fixed_momP = true;
      read(top, "PiPf", fixed_mom.pi_pf);
    }

    if (mom_range.mom_rangeP && fixed_mom.fixed_momP)
      throw string("both fixed and varying momenta ranges specified - only allow one");

    if (! mom_range.mom_rangeP && ! fixed_mom.fixed_momP)
      throw string("neither fixed nor varying momenta ranges specified - must specify one");

  if (paramtop.count("SaveSolutions") == 1)
    {
      XMLReader top(paramtop, "SaveSolutions");

      read(top, "save", saveSoln.saveP);
      read(top, "dir", saveSoln.dir);
      read(top, "pattern", saveSoln.pattern);
    }




 }
  catch(const std::string& e) 
  {
    cerr << "Caught Exception reading XML: " << e << endl;
    exit(1);
  }
}


//! System for linear least squares
/*! System here is rho-pi FF */
class A0PhotonFFSys: public LLSqComponent<EnsemReal>
{
public:
  //! destructor to help with cleanup;
  A0PhotonFFSys(const A0PhotonFFSysParam& param_);

  //! destructor to help with cleanup;
  ~A0PhotonFFSys() {}

  //! Set qsq to use
  void setQsq(const QsqVal& qsq);

  //! Reset time fit interval
  void setTimeFit(int ti, int tf);

  //! Get time fit interval
  ArrayInt getTimeFit() const;

  //! Get time fit ranges
  const Array<ArrayInt>& getFitRange() const;

  //! Given source coordinates, return the corresponding matrix element
  EnsemReal mat(const ArrayInt& ind, int nFF) const;

  //! Given source coordinates, return the corresponding rhs element
  EnsemReal rhs(const ArrayInt& ind) const;

  //! Number of bins
  int nbins() const {return nbin;}

  //! Time extent
  int timeLen() const {return param.lattice.latt_size[param.lattice.decay_dir];}

  //! Max number of FF
  int maxNUnknowns() const {return maxNUnknown;}

  //! Number of FF
  int nUnknowns() const {return nUnknown;}

  //! Set number of FF
  void setNUnknowns(int n);

  //! Number of levels at source
  int nLevelsI() const {return n_levels_i;}

  //! Number of levels at sink
  int nLevelsF() const {return n_levels_f;}

  //! Set number of levels at source
  void setNLevelsI(int n);
  
  //! Set number of levels at source
  void setNLevelsF(int n);
  
  //! Return the dimensions of all the indices
  const ArrayInt& size() const {return sz;}

  //! Energy of initial state at momentum p
  EnsemReal getMassI(int n) const;

  //! Energy of final state at momentum p
  EnsemReal getMassF() const;


  //! Possibly save the solutions
  void save(const EnsemVectorReal& F_array, 
	    const Real& qsq, int t_i, int t_f) const;

protected:
  //! Convenience - make a threept arg
  ThreePtArg threePtArg(const PiPf& pi_pf, int g, int nu) const;

 //! Convenience - make a threept arg with link insertions
  ThreePtArg threePtArg(const PiPf& pi_pf, int g, int nu, int link_dir) const;

  //! Convenience - make a twopt arg
  TwoPtArg twoPtArg(const ArrayInt& p) const;

  //! Energy of rho at momentum p
  EnsemReal energy_i(int level, const ArrayInt& p) const;

  //! Energy of pi at momentum p
  EnsemReal energy_f(const ArrayInt& p) const;

  //! Compute matrix element piece
  EnsemReal matelem_c(int m, const ArrayInt& p_f,  const ArrayInt& p_i, 
		    int mu, int nu) const;
  //! Compute matrix element piece
  EnsemReal matelem_n(int m, const ArrayInt& p_f,  const ArrayInt& p_i, 
		    int mu, int nu) const;
  //! Compute propagation factors
  EnsemReal propfact(int m, const ArrayInt& p_f,  const ArrayInt& p_i, int t) const;

  //! Decompose index
  void decompose(int& m, int& n, int nFF) const;

protected:
  //! Hide default onstructor
  A0PhotonFFSys(); 

private:
  A0PhotonFFSysParam  param;

  mutable Manage3PtFuncBB      threept;
  //  mutable Array<ManageEnergy>  E_photon;
  mutable Array<ManageEnergy>  E_a0;
  // mutable Array<ManageZAmp>    Z_rho;
  mutable Array<ManageZAmp>    Z_a0;

  EnsemReal    Z_V;

  int              nUnknown;
  int              maxNUnknown;
  int              n_levels_i;
  int              n_levels_f;
  int              nbin;
  int              t_i;
  int              t_f;
  Array<PiPf>      pi_pf;
  ArrayInt         sz;
};



// Constructor
A0PhotonFFSys::A0PhotonFFSys(const A0PhotonFFSysParam& param_) : 
  param(param_), 
  threept(new A0Photon3PtFunc(param.threept.pattern), param.threept.cache_file, param.threept.cfg_file)
{
  t_i = param.fit_range[0][0];   //! start of fit range
  t_f = param.fit_range[0][1];   //! end of fit range

  // Read Z_V
  read(param.Z_V_file_name, Z_V);
  nbin = Z_V.size();
  cout << __func__ << ": found nbins=" << nbin << endl;

  // Read energies
  // E_photon.resize(param.E.pattern_photon.size());
  // for(int i=0; i < param.E.pattern_photon.size(); ++i)
  //  E_photon[i].create(new LocalEnergyFunc(param.E.pattern_photon[i], "photon"));

  // Read energies
  E_a0.resize(param.E.pattern_a0.size());
  for(int i=0; i < param.E.pattern_a0.size(); ++i)
    E_a0[i].create(new LocalEnergyFunc(param.E.pattern_a0[i], "a0"));

  // Read amplitudes
  //  Z_rho.resize(param.twopt.pattern_Z_rho.size());
  // for(int i=0; i < param.twopt.pattern_Z_rho.size(); ++i)
  //  Z_rho[i].create(new LocalZFunc(param.twopt.pattern_Z_rho[i], "rho"), param.twopt.avg_2pt_func);

  // Read amplitudes
  Z_a0.resize(param.twopt.pattern_Z_a0.size());
  for(int i=0; i < param.twopt.pattern_Z_a0.size(); ++i)
    Z_a0[i].create(new LocalZFunc(param.twopt.pattern_Z_a0[i], "a0"), param.twopt.avg_2pt_func);

  // Check array sizes - number of exponentials
  if ( Z_a0.size() !=  E_a0.size())
  {
    cerr << __func__ << ": inconsistent array size for E and Z" << endl;
    exit(1);
  }

  // Check ground state energies  
  ArrayInt zero(Nd-1);
  zero = 0;
  if (E_a0[0][zero].size() != nbins())
  {
    cerr << __func__ << ": energies inconsistent size with Z_V" << endl;
    exit(1);
  }

  //n_levels_f = E_photon.size();
  n_levels_i = E_a0.size();
//  maxNUnknown = 2 * Nd * E.size()*(E.size() + 1)/2;
  maxNUnknown = 4;   // temporal & spatial version of V
  nUnknown    = maxNUnknown;
  
  cout << __func__ << ": number of energy levels at source = " << nLevelsI() << endl;
  // cout << __func__ << ": number of energy levels at sink   = " << nLevelsF() << endl;
  cout << __func__ << ": maximum number of unknowns = " << maxNUnknown << endl;
}


//! Convenience - make a threept arg
ThreePtArg A0PhotonFFSys::threePtArg(const PiPf& pi_pf, int g, int nu) const
{
  ThreePtArg  arg;
  arg.pi_pf  = pi_pf;
  arg.g      = g;
  arg.snk.resize(1);
  arg.snk[0] = nu;

//  std::cout << "3pt[]:  g= " << arg.g 
//	    << "   p_f= " << arg.pi_pf.p_f 
//	    << "   p_i= " << arg.pi_pf.p_i
//	    << "   nu= "  << nu
//	    << std::endl;

  return arg;
}

//! Convenience - make a threept arg with link insertions
ThreePtArg A0PhotonFFSys::threePtArg(const PiPf& pi_pf, int g, int nu, int link_dir) const
{
  ThreePtArg  arg;
  arg.pi_pf  = pi_pf;
  arg.g      = g;
  arg.snk.resize(1);
  arg.snk[0] = nu;
  arg.links.resize(1);
  arg.links[0] = link_dir;

  //   std::cout << "3pt[]:  g= " << arg.g 
  //	    << "   p_f= " << arg.pi_pf.p_f 
  //	    << "   p_i= " << arg.pi_pf.p_i
  //    << "   nu= "  << nu
  //     << " link_dir= " << link_dir
  //    << std::endl;

  return arg;
}


  //! Convenience - make a twopt arg
TwoPtArg A0PhotonFFSys::twoPtArg(const ArrayInt& p) const
{
  TwoPtArg  arg;
  arg.p = p;
  return arg;
}



//! Set number of FF
void A0PhotonFFSys::setNUnknowns(int n) 
{
  if (n > maxNUnknowns())
  {
    cerr << __func__ << ": too many unknowns\n";
    exit(1);
  }
  nUnknown = n;
}


// Reset Qsq to be used
void A0PhotonFFSys::setQsq(const QsqVal& qsq)
{
  pi_pf.resize(qsq.pi_pf.size());

  int n = 0;
  for(list<PiPf>::const_iterator nqq=qsq.pi_pf.begin(); 
      nqq != qsq.pi_pf.end(); 
      ++nqq,++n)
  {
    pi_pf[n] = *nqq;
  }

  // Index sizes
  sz.resize(4);
  sz[0] = t_f - t_i + 1;   // fitting time range
  sz[1] = Nd;              // directions
  sz[2] = Nd;            // number of spatial directions
  sz[3] = pi_pf.size();    // number of momenta pairs
}


//! Reset time fit interval
void A0PhotonFFSys::setTimeFit(int ti, int tf) 
{
  t_i = ti; t_f = tf;

  sz[0] = t_f - t_i + 1;   // reset time range
}


//! Get time fit interval
ArrayInt A0PhotonFFSys::getTimeFit() const
{
  ArrayInt tt(2);
  tt[0] = t_i;
  tt[1] = t_f;

  return tt;
}

//! Get time fit interval
const Array<ArrayInt>& A0PhotonFFSys::getFitRange() const
{
  return param.fit_range;
}


//! Set number of levels at source
void A0PhotonFFSys::setNLevelsI(int n)
{
  if (n > E_a0.size())
  {
    cerr << __func__ << ": number of levels too high\n";
    exit(1);
  }
  n_levels_i = n;
}
  
//! Set number of levels at source
//void A0PhotonFFSys::setNLevelsF(int n)
//{
//  if (n > E_photon.size())
//  {
//    cerr << __func__ << ": number of levels too high\n";
//    exit(1);
//  }
//  n_levels_f = n;
//}
  

//! Compute or return energy
EnsemReal A0PhotonFFSys::energy_i(int level, const ArrayInt& p0) const
{
//  cout << __func__ << endl;

  ArrayInt p = canonicalOrder(p0);
  ArrayInt zero(Nd-1);
  zero = 0;

#if 0
  if (p != zero && ! E_a0[level].exist(p))
    E_a0[level].insert(p, energyDisp(E_a0[level][zero], p, param.lattice));
#endif

  return E_a0[level][p];
}


//! Compute or return energy of photon using right dispn reln
EnsemReal A0PhotonFFSys::energy_f(const ArrayInt& p) const
{
//  cout << __func__ << endl;
     const double xi        =  param.lattice.aniso.xi;
     const double xi_sq     =  xi * xi;

     Array<double> nnf = contMom(p,param.lattice.latt_size[0]);
       Array<Real> nf(nnf.size());
      for(int i=0; i < nf.size(); ++i)
        nf[i] = nnf[i];
     const Real test = ( nf[0]*nf[0] + nf[1]*nf[1] + nf[2]*nf[2] )/Real(xi_sq) - Real(param.E.photon_Qsq);

    EnsemReal e ;
    e.resize(nbins());

        if (toBool(  Real(test) >= Real(0.0)  ))
        {
	 e =  photonenergyDisp(getMassF(), p, param.lattice);
        }
		   else
        {
	      e = Real(0.0);
        }


	 return e ;
}


//! Energy of initial state at momentum p
EnsemReal A0PhotonFFSys::getMassI(int n) const
{
  ArrayInt zero(Nd-1);
  zero = 0;
  return E_a0[n][zero];
}


//! Energy of final state at momentum p
EnsemReal A0PhotonFFSys::getMassF() const
{
  //actually the Qsq value for the LSZ photon
  EnsemReal Qsq;
  //    cout << "\nbins in mass_f  " << nbins() << endl;
  Qsq.resize(nbins());
  Qsq = param.E.photon_Qsq;
  return Qsq;

}


//! Compute matrix element piece
EnsemReal A0PhotonFFSys::matelem_c(int m, const ArrayInt& p_f,  const ArrayInt& p_i, 
				int mu, int nu) const

  // current conserving kinematic factor

{
  if (nu > Nd-1)
  {
    cerr << __func__ << ": currently only support mu=0..3, nu<3\n";
    exit(1);
  }

  const double xi  = param.lattice.aniso.xi;
  const double csq = param.lattice.aniso.c_sq;
  ArrayDouble pp_i = sqrt(csq) * contMom(p_i,param.lattice.latt_size[0]) / xi;
  ArrayDouble pp_f = sqrt(csq) * contMom(p_f,param.lattice.latt_size[0]) / xi;

  ArrayInt zero(Nd-1);
  zero = 0;

  const EnsemReal  m_f = getMassF();
  const EnsemReal  m_i = getMassI(m);
  const EnsemReal  E_f = energy_f(p_f);
  const EnsemReal  E_i = energy_i(m,p_i);
  EnsemReal  q2   = minkQsq(E_i, pp_i, E_f, pp_f);
  EnsemReal omega = Real(0.25)*(pow(m_f - m_i,2) - q2)*(pow(m_f + m_i,2) - q2);
  EnsemReal    dt = dot(E_f,pp_f,E_i,pp_i);

  // Transition matrix element
  EnsemReal e;

  if (mu == Nd-1)
  {
      // Temporal
    e.resize(nbins());  
    if (nu != Nd-1)
      {
	e = - (Real(sqrt(xi*xi*xi))/m_i) * E_f * ( Real(pp_i[nu]) - Real(pp_f[nu]) )  ;
      }
    else
      {
	 e = ( Real(sqrt(xi*xi*xi)) / m_i ) * (  (dt + m_f)  -  E_f * ( E_i - E_f ) )    ;
      }

  }
  else if (mu < Nd-1)
  {
    Real dd = Real((mu == nu) ? 1.0 : 0.0);
    // Spatial
 
    e.resize(nbins()); 
    //NB m_f is the Qsq of the LSZ photon
 
    if (nu != Nd-1)
      {
	e =  ( Real(sqrt(xi*xi*xi)) / m_i ) * ( - dd * (dt + m_f) - Real(pp_f[mu]) * (Real(pp_i[nu]) - Real(pp_f[nu]) ) );
       }
    else
      {
	e =  ( Real(sqrt(xi*xi*xi)) / m_i ) * ( - Real(pp_f[mu])*(E_i - E_f)   );
      }
  }
  else 
  {
    cerr << __func__ << ": unexpected direction\n";
    exit(1);
  }

  return e;
}


//! Compute matrix element piece
EnsemReal A0PhotonFFSys::matelem_n(int m, const ArrayInt& p_f,  const ArrayInt& p_i, 
				int mu, int nu) const

{
  if (nu > Nd-1)
  {
    cerr << __func__ << ": currently only support mu=0..3, nu<3\n";
    exit(1);
  }

  const double xi  = param.lattice.aniso.xi;
  const double csq = param.lattice.aniso.c_sq;
  ArrayDouble pp_i = sqrt(csq) * contMom(p_i,param.lattice.latt_size[0]) / xi;
  ArrayDouble pp_f = sqrt(csq) * contMom(p_f,param.lattice.latt_size[0]) / xi;

  ArrayInt zero(Nd-1);
  zero = 0;

  const EnsemReal  m_f = getMassF();
  const EnsemReal  m_i = getMassI(m);
  const EnsemReal  E_f = energy_f(p_f);
  const EnsemReal  E_i = energy_i(m,p_i);
  EnsemReal  q2   = minkQsq(E_i, pp_i, E_f, pp_f);
  EnsemReal omega = Real(0.25)*(pow(m_f - m_i,2) - q2)*(pow(m_f + m_i,2) - q2);
  EnsemReal    dt = dot(E_f,pp_f,E_i,pp_i);

  // Transition matrix element
  EnsemReal e;

  if (mu == Nd-1)
  {
      // Temporal
    e.resize(nbins());
    if (nu != Nd-1)
      {
	e = (Real(sqrt(xi*xi*xi))/m_i) * E_f * ( Real(pp_i[nu]) - Real(pp_f[nu]) )    ;
      }
    else 
      {
	e = ( Real(sqrt(xi*xi*xi)) / m_i ) * (  (dt + m_f)  +  E_f * ( E_i - E_f ) );
      }
  }
  else if (mu < Nd-1)
  {
    Real dd = Real((mu == nu) ? 1.0 : 0.0);
    // Spatial
 
    e.resize(nbins()); 
    //NB m_f is the Qsq of the LSZ photon
 
    if (nu != Nd-1)
      {
	e =  ( Real(sqrt(xi*xi*xi)) / m_i ) * ( - dd * (dt + m_f) + Real(pp_f[mu]) * (Real(pp_i[nu]) - Real(pp_f[nu])   ));
       }
    else
      {
	e =  ( Real(sqrt(xi*xi*xi)) / m_i ) * ( Real(pp_f[mu])*(E_i - E_f)   );
      }
  }

  else 
  {
    cerr << __func__ << ": unexpected direction\n";
    exit(1);
  }

  return e;
}


//! Compute propagation factors
EnsemReal A0PhotonFFSys::propfact(int m, const ArrayInt& p_f, const ArrayInt& p_i, int t) const
{
  int T_i = param.t_source;
  int T_f = param.t_sink;
  int T   = param.lattice.latt_size[param.lattice.decay_dir];

  EnsemReal g, g_i, g_f;
  g.resize(nbins());
  g_i.resize(nbins());
  g_i = exp(-Real(t - T_i)*energy_i(m,p_i));
  g_f = exp( energy_f(p_f)*Real(t));
  g  = Z_a0[m][twoPtArg(p_i)] * g_i * g_f;
  g /= Real(2) * energy_i(m,p_i);

  return g;
}


//! Decompose index
/*! A really dumb way to find coordinates within a lower triangular matrix */
void A0PhotonFFSys::decompose(int& m, int& n, int nFF) const
{
  // Just run through the goofy thing
  for(n=0; n < nLevelsI(); ++n)
    for(m=0; m <= min(nLevelsF()-1,n); ++m)
    {
      int tt = n*(n+1)/2 + m;
      if (nFF == tt)
	return;
    }

  cerr << __func__ << ": failed to decompose index: nFF=" << nFF << endl;
  exit(1);
}

// Given source coordinates, return the corresponding matrix element
EnsemReal A0PhotonFFSys::mat(const ArrayInt& ind, int nFF) const
{
//  cout << __func__ << endl;

  // NOTE: encoding here:   nFF has indices like  int[numTemporalFF][2]
  // So, if m and n encode the energy level, then one indexes like
  //
  // (nFF & 1) == 0 spatial V
  // (nFF & 3) == 1 temporal V
  // nFF = <spatial/temporal V> + 2*( n*n + m)
  //
  int t  = ind[0] + t_i;
  int mu = ind[1];
  int nu = ind[2];
  int p  = ind[3];

  ArrayInt p_i = pi_pf[p].p_i;
  ArrayInt p_f = pi_pf[p].p_f;

  int m = 0;
  //decompose(m, n, nFF >> 1); // Determine sink,source energy levels

  EnsemReal f;
  f.resize(nbins());
  f = 0;

  // NOTE: should eventually support cross term in transition FF with additional propfacts

  // spatial C
  if ((nFF) == 0 && (mu < Nd-1))
  {
    f = matelem_c(m, p_f, p_i, mu, nu) * propfact(m, p_f,  p_i, t);
  }

  // temporal C
  if ((nFF) == 1 && (mu == Nd-1))
  {
    f = matelem_c(m, p_f, p_i, mu, nu) * propfact(m, p_f,  p_i, t);
  }
  // spatial N
    if ((nFF) == 2 && (mu < Nd-1))
   {
    f = matelem_n(m, p_f, p_i, mu, nu) * propfact(m, p_f,  p_i, t);
   }

  // temporal N
   if ((nFF) == 3 && (mu == Nd-1))
   {
    f = matelem_n(m, p_f, p_i, mu, nu) * propfact(m, p_f,  p_i, t);
   }
#if 0
  char ss[1024];
  sprintf(ss,"mat_f%d_mu%d_nu%d_m%d_pf%d%d%d_n%d_pi%d%d%d", nFF,mu,nu,m,p_f[0],p_f[1],p_f[2],n,p_i[0],p_i[1],p_i[2]);
  write(string(ss),f);
#endif

  // Final kinetic factor
  f /= Z_V*Z_V;

//  cout << "exiting mat" << endl;

  return f;
}


// Given source coordinates, return the corresponding rhs element
EnsemReal A0PhotonFFSys::rhs(const ArrayInt& ind) const
{
//  cout << __func__ << endl;

  int t  = ind[0] + t_i;
  int mu = ind[1];
  int nu = ind[2];
  int p  = ind[3];

  if (nu < 0 || nu >= Nd-1)
  {
    //cerr << __func__ << ": nu out of bounds\n";
    // exit(1);
  }

  //point-like vector current
  // EnsemVectorComplexF& three = threept[threePtArg(pi_pf[p], 1 << mu, nu)];
  EnsemReal thr;
  EnsemVectorComplexF three;

  //we can construct the conserved current here out of the one-link bb's
  if ( mu == Nd-1) //temporally directed current
   {
     EnsemVectorComplexF& three_a = threept[threePtArg(pi_pf[p], 0, nu, 3)];
      EnsemVectorComplexF& three_b = threept[threePtArg(pi_pf[p], 8, nu, 3)];
     EnsemVectorComplexF& three_c = threept[threePtArg(pi_pf[p], 0, nu, 7)];
     EnsemVectorComplexF& three_d = threept[threePtArg(pi_pf[p], 8, nu, 7)];

     three = Real(0.5) * ( shift( three_a + three_b , 1 ) - ( three_c - three_d ) );

   }
  else
   {
     EnsemVectorComplexF& three_a = threept[threePtArg(pi_pf[p], 0, nu, mu )];
     EnsemVectorComplexF& three_b = threept[threePtArg(pi_pf[p], 1 << mu, nu, mu )];
     EnsemVectorComplexF& three_c = threept[threePtArg(pi_pf[p], 0, nu, mu + 4 )];
     EnsemVectorComplexF& three_d = threept[threePtArg(pi_pf[p], 1 << mu, nu, mu + 4 )];

     ArrayInt p_i = pi_pf[p].p_i;
     ArrayInt p_f = pi_pf[p].p_f;
     const double xi  = param.lattice.aniso.xi;
     const double csq = param.lattice.aniso.c_sq;
     ArrayDouble pp_i = sqrt(csq) * contMom(p_i,param.lattice.latt_size[0]) / xi;
     ArrayDouble pp_f = sqrt(csq) * contMom(p_f,param.lattice.latt_size[0]) / xi; 
 
     Complex phase = cmplx(  cos( Real(pp_f[mu]) - Real(pp_i[mu])), -sin( Real(pp_f[mu]) - Real(pp_i[mu])) );
      three = Real(0.5) * ( phase*(three_a + three_b) - ( three_c - three_d ) );    
    }

  int sg = ((nu == 3) && (mu == 3) )? -1 : 1;

  // Extract appropriate part of 3pt. Fold in rho_y sign above.
  // the minus sign corrects for the use of the euclidean antisymmetric four-tensor
  if (  ((mu == Nd-1) && (nu != Nd-1)) || ((mu != Nd-1) && (nu == Nd-1))  )
    thr =  Real(sg)*imag(peekObs(three,t)) ;
  else 
    thr =  Real(sg)*real(peekObs(three,t));

#if 0
  ArrayInt p_i = pi_pf[p].p_i;
  ArrayInt p_f = pi_pf[p].p_f;

  char ss[1024];
  sprintf(ss,"rhs_mu%d_nu%d_pf%d%d%d_pi%d%d%d", mu,nu,p_f[0],p_f[1],p_f[2],p_i[0],p_i[1],p_i[2]);
  write(string(ss),thr);
#endif

//  cout << "exiting rhs" << endl;

  return thr;
}


//! Possibly save the solutions
void A0PhotonFFSys::save(const EnsemVectorReal& F_array, 
			const Real& qsq, int t_i, int t_f) const
{
 if (param.saveSoln.saveP)
 {
   char dir[1024];
   sprintf(dir, param.saveSoln.dir.c_str(), 
	    param.mom_range.sink_mom[0][0], 
	    param.mom_range.sink_mom[0][1], 
	    param.mom_range.sink_mom[0][2], 
	    toFloat(qsq));

   cout << "mkdir: " << dir << endl;
   mkdir(dir, 0755);

    for(int n=0; n < nUnknowns(); ++n)
    {
      char ss[1024];
      EnsemReal F = peekObs(F_array, n);

      sprintf(ss, param.saveSoln.pattern.c_str(),
	      n, t_i, t_f);

      string filename = string(dir) + '/' + string(ss);
      
      write(filename, F);
    }
  }
}



int main(int argc, char *argv[])
{
  if (argc != 3)
  {
    cerr << "Usage:"
	 << argv[0] << ":  <input xml file>   <output xml file>" << endl;
    exit(1);
  }

  XMLReader xml_in(argv[1]);
  A0PhotonFFSysParam linParam(xml_in, "/A0PhotonFit");

  XMLFileWriter xml_out(argv[2]);
  push(xml_out, "A0PhotonFF");

  cout << "construct linear system" << endl;
  A0PhotonFFSys sys(linParam);

  cout << "sys.nbins=" << sys.nbins() << endl;

  // The momenta list
  std::list<QsqVal> qsq_val;

  // Use some method to construct Qsq list
  if (linParam.mom_range.mom_rangeP)
  {
    cout << "Construct qsq list over momenta range" << endl;
    QsqParam qsqParam;
    qsqParam.nq_norm_sq_max = linParam.mom_range.mom2_max;
    qsqParam.lattice        = linParam.lattice;
    qsqParam.mass_i         = toDouble(mean(sys.getMassI(0)));
    qsqParam.mass_f         = toDouble(mean(sys.getMassF())); //now the LSZmass of the photon
  
    for(int i=0; i < linParam.mom_range.sink_mom.size(); ++i)
      qsqParam.nf_list.push_back(linParam.mom_range.sink_mom[i]);

    // Construct list
    qsq_val = constructQsq(qsqParam);

#if 1
    // Cannot use 0 momentum
    //   {
    //     ArrayInt zero(Nd-1);
    //    zero = 0;
    //    if (qsq_val.front().pi_pf.front().p_f == zero &&
    //	  qsq_val.front().pi_pf.front().p_i == zero)
    //   {
    //	cout << "First momenta has p_f=p_i=zero - drop since it will cause problems\n" << endl;
    //	qsq_val.pop_front();
    //    }
    //  }
#endif

    // Print out list
    {
      std::ofstream f(linParam.qsq_file_name.c_str());
      f << qsq_val;
      f.close();
    }
  }
  else if (linParam.fixed_mom.fixed_momP)
  {
    cout << "Use fixed momentum list" << endl;

    ArrayInt p_f = linParam.fixed_mom.pi_pf[0].p_f;
    ArrayInt p_i = linParam.fixed_mom.pi_pf[0].p_i;
    const double mass_f = toDouble(mean(sys.getMassF()));
    const double mass_i = toDouble(mean(sys.getMassI(0)));

    QsqVal  qsq;
    qsq.Qsq_c = -qsqContDispPhys(mass_f, p_f, mass_i, p_i, linParam.lattice);
    qsq.Qsq_l = -qsqLattDispPhys(mass_f, p_f, mass_i, p_i, linParam.lattice);
    qsq.qsq   = norm2(p_f - p_i);

    for(int i=0; i < linParam.fixed_mom.pi_pf.size(); ++i)
      qsq.pi_pf.push_back(linParam.fixed_mom.pi_pf[i]);

    qsq_val.push_back(qsq);
  }
  else
  {
    cerr << argv[0] << ": invalid momentum specification\n";
    exit(1);
  }

  push(xml_out, "VaryUnknowns");
  int i = sys.nUnknowns();
  {
    push(xml_out, "elem");
    sys.setNUnknowns(i);
    cout << "\nNumber of unknowns set to  " << sys.nUnknowns() << endl;
    cerr << "\nNumber of unknowns set to  " << sys.nUnknowns() << endl;
    ff_driver(xml_out, sys, qsq_val);
    pop(xml_out);
  }
  pop(xml_out);

  pop(xml_out);  // A0PhotonFF
  xml_out.close();
}
