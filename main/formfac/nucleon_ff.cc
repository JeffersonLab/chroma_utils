// $Id: nucleon_ff.cc,v 2.0 2008/12/05 04:43:48 edwards Exp $
// Nucleon form-factor code - for LHPC analysis. No anisotropy

#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "formfac/formfac.h"
#include <list>
#include <iostream>

using namespace std;

namespace FF
{
  //! Driver
  void ff_driver(XMLWriter& xml_out, LLSqComponent<EnsemComplex>& sys, const std::list<QsqVal>& qsq_val);

  //-----------------------------------------------------------------------------------

  //! Make a state
  class Local3PtFunc : public State3PtFunc
  {
  public:
    //! Constructor
    Local3PtFunc(const std::string& pattern_) : pattern(pattern_) {}

    //! Construct 3pt names
    std::string operator()(int cfg, const ThreePtArg& arg) const
      {
	ArrayInt p_i = arg.pi_pf.p_i;
	ArrayInt p_f = arg.pi_pf.p_f;
	ArrayInt q   = p_f - p_i;      // note sign convention
  
//	cout << __func__ << ": pattern=" << pattern.c_str() << std::endl;

	char qxsign = (q[0] < 0 ) ? '-' : '+';
	char qysign = (q[1] < 0 ) ? '-' : '+';
	char qzsign = (q[2] < 0 ) ? '-' : '+';

	char pfxsign = (p_f[0] < 0 ) ? '-' : '+';
	char pfysign = (p_f[1] < 0 ) ? '-' : '+';
	char pfzsign = (p_f[2] < 0 ) ? '-' : '+';

	if (arg.snk[0] < 0 || arg.snk[0] >= 2)
	{
	  std::cerr << __func__ << ": snk out of bounds\n";
	  exit(1);
	}
	char ud[2] = {'U', 'D'};

	char filen[1024];
	sprintf(filen, pattern.c_str(), 
		cfg,
		ud[arg.snk[0]],
		pfzsign, abs(p_f[2]), pfysign, abs(p_f[1]), pfxsign, abs(p_f[0]),
		qzsign, abs(q[2]), qysign, abs(q[1]), qxsign, abs(q[0]));

	cout << "3pt filen=" << filen << endl;

	return std::string(filen);
      }

  private:
    std::string pattern;
  };

  
  //! State function
  class Local2PtFunc : public State2PtFunc
  {
  public:
    //! Constructor
    Local2PtFunc(const std::string& pattern_, const std::string& particle_) :
      pattern(pattern_), particle(particle_) {}

    //! Construct 2pt names
    std::string operator()(const TwoPtArg& arg) const
      {
	ArrayInt zero(Nd-1);
	zero = 0;
	std::ostringstream pion;

	pion << particle;

	if (arg.p != zero)
	  pion << "_px" << arg.p[0] << "_py" << arg.p[1] << "_pz" << arg.p[2];

	char filen[1024];
	sprintf(filen, pattern.c_str(), pion.str().c_str());

	cout << "2pt filen=" << filen << endl;

	return std::string(filen);
      }

  private:
    std::string pattern;
    std::string particle;
  };


  //! State function
  class LocalZFunc : public State2PtFunc
  {
  public:
    //! Constructor
    LocalZFunc(const std::string& pattern_, const std::string& particle_) :
      pattern(pattern_), particle(particle_) {}

    //! Construct Z names
    std::string operator()(const TwoPtArg& arg) const
      {
	ArrayInt zero(Nd-1);
	zero = 0;
	std::ostringstream pion;

	pion << particle;

	if (arg.p != zero)
	  pion << "_px" << arg.p[0] << "_py" << arg.p[1] << "_pz" << arg.p[2];

	char filen[1024];
	sprintf(filen, pattern.c_str(), pion.str().c_str());

	return std::string(filen);
      }

  private:
    std::string pattern;
    std::string particle;
  };


  //! State function
  class LocalEnergyFunc : public StateEnergyFunc
  {
  public:
    //! Constructor
    LocalEnergyFunc(const std::string& pattern_, const std::string& particle_) : 
      pattern(pattern_), particle(particle_) {}

    // Construct energy names
    std::string operator()(const ArrayInt& mom) const
      {
	ArrayInt zero(Nd-1);
	zero = 0;
	std::ostringstream pion;

	pion << particle;

	if (mom != zero)
	  pion << "_px" << mom[0] << "_py" << mom[1] << "_pz" << mom[2];

	char filen[1024];
	sprintf(filen, pattern.c_str(), pion.str().c_str());

	return std::string(filen);
      }

  private:
    std::string pattern;
    std::string particle;
  };


} // namespace FF


//using namespace ADATXML;
using namespace ENSEM;
using namespace FF;



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
struct NucleonFFSysParam
{
  NucleonFFSysParam();
  NucleonFFSysParam(XMLReader& xml_in, const std::string& path);

  struct Manage3PtParam
  {
    std::string  pattern;
    std::string  cache_file;
    std::string  cfg_file;
    int          max_map_megabytes;
  } threept;

  struct Manage2PtParam
  {
    std::string   pattern_SP;
    std::string   pattern_SS;
    bool          avg_2pt_func;
  } twopt;

  struct ManageEnergyParam
  {
    std::string   pattern;
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
};

// Nucleon sys params
NucleonFFSysParam::NucleonFFSysParam(XMLReader& xml_in, const std::string& path) 
{
  try 
  {
    XMLReader paramtop(xml_in, path);

    {
      XMLReader top(paramtop, "ThreePt");

      read(top, "pattern", threept.pattern);
      read(top, "cache_file", threept.cache_file);
      read(top, "cfg_file", threept.cfg_file);
      read(top, "max_map_megabytes", threept.max_map_megabytes);
    }

    {
      XMLReader top(paramtop, "TwoPt");

      read(top, "pattern_SP", twopt.pattern_SP);
      read(top, "pattern_SS", twopt.pattern_SS);
      read(top, "avg_2pt_func", twopt.avg_2pt_func);
    }

    {
      XMLReader top(paramtop, "Energy");

      read(top, "pattern", E.pattern);
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
  }
  catch(const std::string& e) 
  {
    cerr << "Caught Exception reading XML: " << e << endl;
    exit(1);
  }
}


//! System for linear least squares
/*! System here is nucleon FF */
class NucleonFFSys: public LLSqComponent<EnsemComplex>
{
public:
  //! destructor to help with cleanup;
  NucleonFFSys(const NucleonFFSysParam& param_);

  //! destructor to help with cleanup;
  ~NucleonFFSys() {}

  //! Set qsq to use
  void setQsq(const QsqVal& qsq);

  //! Reset time fit interval
  void setTimeFit(int ti, int tf) {t_i = ti; t_f = tf;}

  //! Get time fit interval
  ArrayInt getTimeFit() const;

  //! Get time fit ranges
  const Array<ArrayInt>& getFitRange() const;

  //! Given source coordinates, return the corresponding matrix element
  EnsemComplex mat(const ArrayInt& ind, int nFF) const;

  //! Given source coordinates, return the corresponding rhs element
  EnsemComplex rhs(const ArrayInt& ind) const;

  //! Number of bins
  int nbins() const {return nbin;}

  //! Time extent
  int timeLen() const {return param.lattice.latt_size[param.lattice.decay_dir];}

  //! Number of FF
  /*! F_1^p -> 0; F_2^p -> 1, F_1^n -> 2, F_2^n -> 3 */
  int nUnknowns() const {return 4;}

  //! Return the dimensions of all the indices
  const ArrayInt& size() const {return sz;}

  //! Energy of initial state at momentum p
  EnsemReal getMassI() const;

  //! Energy of final state at momentum p
  EnsemReal getMassF() const;

protected:
  //! Convenience - make a threept arg
  ThreePtArg threePtArg(const PiPf& pi_pf, int g, int ud) const;

  //! Convenience - make a twopt arg
  TwoPtArg twoPtArg(const ArrayInt& p) const;

  //! Energy of initial state at momentum p
  EnsemReal energyI(const ArrayInt& p) const;

  //! Energy of final state at momentum p
  EnsemReal energyF(const ArrayInt& p) const;

  //! Return spin projector in use
  SpinMatrix projectorT() const;
  
  //! Compute a fermion term at source (same as sink)
  EnsemSpinMatrix fermI(const ArrayInt& p_i) const;

  //! Compute a fermion term at sink (same as source)
  EnsemSpinMatrix fermF(const ArrayInt& p_f) const;

  //! Compute matrix element piece
  EnsemComplex matelem_f1(const ArrayInt& p_f, const ArrayInt& p_i, 
			  int mu) const;
  
  //! Compute matrix element piece
  EnsemComplex matelem_f2(const ArrayInt& p_f, const ArrayInt& p_i, 
			  int mu) const;

protected:
  //! Hide default onstructor
  NucleonFFSys(); 

private:
  NucleonFFSysParam  param;

  mutable ManageEnergyMap      E;
  mutable Manage2PtFuncMap     twopt_SP;
  mutable Manage2PtFuncMap     twopt_SS;
  mutable Manage3PtFuncCacheBB threept;

  EnsemReal    Z_V;

  int          t_i, t_f;
  int          nbin;
  Array<PiPf>  pi_pf;
  ArrayInt     sz;
};



// Constructor
NucleonFFSys::NucleonFFSys(const NucleonFFSysParam& param_) : 
  param(param_), E(new LocalEnergyFunc(param.E.pattern, "proton"), param.lattice), 
  twopt_SP(new Local2PtFunc(param.twopt.pattern_SP, "proton"), param.twopt.avg_2pt_func, param.lattice), 
  twopt_SS(new Local2PtFunc(param.twopt.pattern_SS, "proton"), param.twopt.avg_2pt_func, param.lattice), 
  threept(new Local3PtFunc(param.threept.pattern), 
	  param.threept.cache_file, param.threept.cfg_file, param.lattice, param.threept.max_map_megabytes)
{
  t_i = param.fit_range[0][0];   //! start of fit range
  t_f = param.fit_range[0][1];   //! end of fit range

  // Read ground state energies  
  ArrayInt zero(Nd-1);
  zero = 0;
  nbin = E[zero].size();
  cout << __func__ << ": found nbins=" << nbin << endl;

  // Read Z_V
  read(param.Z_V_file_name, Z_V);
  if (Z_V.size() != nbins())
  {
    cerr << __func__ << ": Z_V inconsistent size with ground state energy" << endl;
    exit(1);
  }

  // NOTE: compensate for technically incorrect Z_V here. I should take into
  // account the BC - the charge on the other side of the sink.
//  Z_V /= Real(0.5);
}


//! Convenience - make a threept arg
ThreePtArg NucleonFFSys::threePtArg(const PiPf& pi_pf, int g, int ud) const
{
  ThreePtArg  arg;
  arg.pi_pf = pi_pf;
  arg.g     = g;
  arg.snk.resize(1);
  arg.snk[0] = ud;

  //  cout << __func__ << "[g=" << g << ", ud=" << ud 
  //       << ", p_f="  << pi_pf.p_f << " , p_i=" << pi_pf.p_i << "]\n";

  return arg;
}

  //! Convenience - make a twopt arg
TwoPtArg NucleonFFSys::twoPtArg(const ArrayInt& p) const
{
  TwoPtArg  arg;
  arg.p = p;

  //  cout << __func__ << "[p="  << p << "]\n";

  return arg;
}



// Reset Qsq to be used
void NucleonFFSys::setQsq(const QsqVal& qsq)
{
  pi_pf.resize(qsq.pi_pf.size());

  int n = 0;
  for(list<PiPf>::const_iterator nqq=qsq.pi_pf.begin(); 
      nqq != qsq.pi_pf.end(); 
      ++nqq,++n)
  {
    pi_pf[n] = *nqq;
  }

  // Clear out any unused BB's - they'll probably not be used at a new qsq
  threept.eraseUnused();

  // Index sizes
  sz.resize(4);
  sz[0] = t_f - t_i + 1;   // fitting time range
  sz[1] = Nd;              // directions
 // sz[1] = Nd-1;              // space directions
 //  sz[1] = 1;              // directions, only Nd-1
  sz[2] = pi_pf.size();    // number of momenta pairs
  sz[3] = 2;               // 0=proton, 1=neutron
}


//! Get time fit interval
ArrayInt NucleonFFSys::getTimeFit() const
{
  ArrayInt tt(2);
  tt[0] = t_i;
  tt[1] = t_f;

  return tt;
}

//! Get time fit interval
const Array<ArrayInt>& NucleonFFSys::getFitRange() const
{
  return param.fit_range;
}


//! Energy of initial state at momentum p
EnsemReal NucleonFFSys::energyI(const ArrayInt& p0) const
{
//  cout << __func__ << endl;

  ArrayInt p = canonicalOrder(p0);
  ArrayInt zero(Nd-1);
  zero = 0;

#if 1
  if (p != zero && ! E.exist(p))
    E.insert(p, energyDisp(E[zero], p, param.lattice));
#endif

  return E[p];
}


//! Energy of final state at momentum p
EnsemReal NucleonFFSys::energyF(const ArrayInt& p0) const
{
  return energyI(p0);
}


//! Energy of initial state at momentum p
EnsemReal NucleonFFSys::getMassI() const
{
  ArrayInt zero(Nd-1);
  zero = 0;
  return E[zero];
}


//! Energy of final state at momentum p
EnsemReal NucleonFFSys::getMassF() const
{
  ArrayInt zero(Nd-1);
  zero = 0;
  return E[zero];
}


//! Compute a fermion term at source (same as sink)
EnsemSpinMatrix NucleonFFSys::fermI(const ArrayInt& p_i) const
{
  const double xi  = param.lattice.aniso.xi;
  ArrayDouble pp_i = contMom(p_i,param.lattice.latt_size[0]) / xi;

  SpinMatrix one = 1.0;
  const EnsemReal& E_i = energyI(p_i);
  EnsemSpinMatrix f = getMassI();

  for(int mu=0; mu < pp_i.size(); ++mu)
    f -= timesI(Real(pp_i[mu]))*(Gamma(1 << mu) * one);

  f += E_i * (Gamma(8) * one);

  return f;
}


//! Compute a fermion term at sink (same as source)
EnsemSpinMatrix NucleonFFSys::fermF(const ArrayInt& p_f) const
{
  return fermI(p_f);
}


//! Yep, sigma_{\mu\nu}
SpinMatrix sigma_munu(int mu, int nu)
{
  SpinMatrix one = 1.0;
  SpinMatrix mat = Real(0.5)* timesI(Gamma(1 << mu)*(Gamma(1 << nu)*one) - 
				     Gamma(1 << nu)*(Gamma(1 << mu)*one));
  return mat;
}



// Projector in use
SpinMatrix NucleonFFSys::projectorT() const
{
  // This is Tmixed() from  chroma/lib/meas/hadron/barspinmats_w.cc

  //! T = (1 + \Sigma_3)*(1 + gamma_4) / 2   = (1 + Gamma(8) - i G(3) - i G(11)) / 2
  SpinMatrix one = 1.0;
  return SpinMatrix(Real(0.5)*(one + Gamma(8)*one + timesMinusI(Gamma(3)*one + Gamma(11)*one)));
}


//! Compute matrix element piece
EnsemComplex NucleonFFSys::matelem_f1(const ArrayInt& p_f, const ArrayInt& p_i, 
				      int mu) const
{
  EnsemComplex f = trace(projectorT() * fermF(p_f) * Gamma(1 << mu) * fermI(p_i));
  return f;
}


//! Compute matrix element piece
EnsemComplex NucleonFFSys::matelem_f2(const ArrayInt& p_f, const ArrayInt& p_i, 
				   int mu) const
{
  // create 4-vector q
  Array<EnsemComplex> q(Nd);
  const EnsemReal& E_f = energyF(p_f);
  const EnsemReal& E_i = energyI(p_i);

  const double xi  = param.lattice.aniso.xi;
  ArrayDouble pp_i = contMom(p_i,param.lattice.latt_size[0]) / xi;
  ArrayDouble pp_f = contMom(p_f,param.lattice.latt_size[0]) / xi;

   q[Nd-1] = timesMinusI(E_f - E_i);
   //  q[Nd-1] = timesI(E_f - E_i);

  for(int j=0; j < Nd-1; ++j)
  {
    q[j].resize(nbins());
    q[j] = Real(pp_f[j]) - Real(pp_i[j]);
  }

  // Compute  \sum sigma_munu * q_nu / (2*m_N);
  EnsemSpinMatrix sum;
  sum.resize(nbins());
  sum = zero;

  sum += sigma_munu(Nd-1,mu) * q[Nd-1];
  for(int nu=0; nu < Nd-1; ++nu)
    sum += sigma_munu(mu,nu) * q[nu];


  sum /= (Real(2.0) * getMassI());

  // Finally compute kinematic factor
  EnsemComplex f = trace(projectorT() * fermF(p_f) * sum * fermI(p_i));
  return f;
}


// Given source coordinates, return the corresponding matrix element
EnsemComplex NucleonFFSys::mat(const ArrayInt& ind, int nFF) const
{
//  cout << __func__ << endl;

  int t    = ind[0] + t_i;
  int mu   = ind[1];
  //  int mu   = Nd-1;
  //  int mu   = 0;
  int p    = ind[2];
  int ud   = ind[3];   // proton=0, neutron=1

  ArrayInt p_i = pi_pf[p].p_i;
  ArrayInt p_f = pi_pf[p].p_f;

  EnsemComplex f;
  f.resize(nbins());
  f = zero;

  // u-quark contribution FF
  if (ud == 0)
  {
    // F_1^u
    if (nFF == 0)
      f = matelem_f1(p_f, p_i, mu);

    // F_2^u
    if (nFF == 1)
      f = matelem_f2(p_f, p_i, mu);
  }

  // d-quark constribution FF
  if (ud == 1)
  {
    // F_1^d
    if (nFF == 2)
      f = matelem_f1(p_f, p_i, mu);

    // F_2^d
    if (nFF == 3)
      f = matelem_f2(p_f, p_i, mu);
  }


#if 0
  {
    char ss[1024];
    sprintf(ss,"mat_f%d_mu%d_ud%d_pf%d%d%d_pi%d%d%d", nFF,mu,ud,p_f[0],p_f[1],p_f[2],p_i[0],p_i[1],p_i[2]);
    write(string(ss),f);
  }
#endif

//  cout << "exiting mat" << endl;

  return f;
}


// Given source coordinates, return the corresponding rhs element
EnsemComplex NucleonFFSys::rhs(const ArrayInt& ind) const
{
//  cout << __func__ << endl;

  int t    = ind[0] + t_i;
  int mu   = ind[1];
  //  int mu   = Nd-1;
//  int mu   = 3;
  int p    = ind[2];
  int ud   = ind[3];   // u=0, d=1

  ArrayInt p_i = pi_pf[p].p_i;
  ArrayInt p_f = pi_pf[p].p_f;

  const EnsemReal E_i = energyI(p_i);
  const EnsemReal m_i = getMassI();

  const EnsemReal E_f = energyF(p_f);
  const EnsemReal m_f = getMassF();

  const EnsemVectorComplexF three_u = threept[threePtArg(pi_pf[p], 1 << mu, 0)];
  const EnsemVectorComplexF three_d = threept[threePtArg(pi_pf[p], 1 << mu, 1)];

  EnsemComplex thr;

  /*  const EnsemVectorComplexF& three = threept[threePtArg(pi_pf[p], 1 << mu, ud)];

  EnsemComplex thr;

  if (mu == Nd-1)
  {
     thr = peekObs(three,t);
 
  }
  else
  {
    thr = peekObs(three,t);    // note i*complex
  }
  */

  if (ud == 0)//u-d
     thr=peekObs(three_u,t)-peekObs(three_d,t);
  if (ud == 1)//u+d
     thr=peekObs(three_u,t)+peekObs(three_d,t);
   

  if (t < param.t_source)
  {
    cerr << __func__ << ": do not currently support peeking 2pts before t_source" << endl;
    exit(1);
  }


  // Ratio
  EnsemComplex f = (thr * peekObs(twopt_SP[twoPtArg(p_f)],t-param.t_source)) /
    (peekObs(twopt_SP[twoPtArg(p_i)],t-param.t_source) * peekObs(twopt_SS[twoPtArg(p_f)],param.t_sink-param.t_source));
  //  f *= Z_V;

  // Correct for 2pt and 3pt energy normalization factors
  f *=  Real(4.0) * E_f * (E_i + m_i);



#if 0
  {
    char ss[1024];
    sprintf(ss,"E_pf%d%d%d", p_f[0],p_f[1],p_f[2]);
    write(string(ss),E_f);

    sprintf(ss,"E_pi%d%d%d", p_i[0],p_i[1],p_i[2]);
    write(string(ss),E_i);

    sprintf(ss,"rhs_mu%d_ud%d_pf%d%d%d_pi%d%d%d", mu,ud,p_f[0],p_f[1],p_f[2],p_i[0],p_i[1],p_i[2]);
    write(string(ss),f);

    sprintf(ss,"twosp_pf_t_mu%d_ud%d_pf%d%d%d_pi%d%d%d", mu,ud,p_f[0],p_f[1],p_f[2],p_i[0],p_i[1],p_i[2]);
    write(string(ss),peekObs(twopt_SP[twoPtArg(p_f)],t-param.t_source));

    sprintf(ss,"twosp_pi_t_mu%d_ud%d_pf%d%d%d_pi%d%d%d", mu,ud,p_f[0],p_f[1],p_f[2],p_i[0],p_i[1],p_i[2]);
    write(string(ss),peekObs(twopt_SP[twoPtArg(p_i)],t-param.t_source));

    sprintf(ss,"twoss_pf_tf_mu%d_ud%d_pf%d%d%d_pi%d%d%d", mu,ud,p_f[0],p_f[1],p_f[2],p_i[0],p_i[1],p_i[2]);
    write(string(ss),peekObs(twopt_SS[twoPtArg(p_f)],param.t_sink-param.t_source));

    sprintf(ss,"thr_mu%d_ud%d_pf%d%d%d_pi%d%d%d", mu,ud,p_f[0],p_f[1],p_f[2],p_i[0],p_i[1],p_i[2]);
    write(string(ss),thr);
  }
#endif

//  cout << "exiting rhs" << endl;

  return f;
}




int main(int argc, char *argv[])
{
  if (argc != 3)
  {
    cerr << "Usage:"
	 << argv[0] << ":  <input xml file>  <output xml file>" << endl;
    exit(1);
  }

  XMLReader xml_in(argv[1]);
  NucleonFFSysParam linParam(xml_in, "/NucleonFF");

  XMLFileWriter xml_out(argv[2]);
  push(xml_out, "NucleonFF");

  cout << "construct linear system" << endl;
  NucleonFFSys sys(linParam);

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
    qsqParam.mass_i         = toDouble(mean(sys.getMassI()));
    qsqParam.mass_f         = toDouble(mean(sys.getMassF()));
  
    for(int i=0; i < linParam.mom_range.sink_mom.size(); ++i)
      qsqParam.nf_list.push_back(linParam.mom_range.sink_mom[i]);

    // Construct list
    qsq_val = constructQsq(qsqParam);

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
    const double mass_i = toDouble(mean(sys.getMassI()));

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

  // Call the linear system driver
  ff_driver(xml_out, sys, qsq_val);

  pop(xml_out);  // NucleonFF
  xml_out.close();
}
