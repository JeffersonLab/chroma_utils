// $Id: pion_2photon.cc,v 2.0 2008/12/05 04:43:49 edwards Exp $
// Pion -> 2photon decay

#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "formfac/formfac.h"
#include "state.h"
#include <list>
#include <iostream>


namespace FF
{
  //! Make a state
  class RhoPion3PtFunc : public State3PtFunc
  {
  public:
    //! Constructor
    RhoPion3PtFunc(const std::string& pattern_) : pattern(pattern_) {}

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

	if (arg.src[0] < 0 || arg.src[0] >= 3)
	{
	  std::cerr << __func__ << ": src out of bounds\n";
	  exit(1);
	}
	char xyz[3] = {'X', 'Y', 'Z'};

	char filen[1024];
	sprintf(filen, pattern.c_str(),
		pfzsign, abs(p_f[2]), pfysign, abs(p_f[1]), pfxsign, abs(p_f[0]),
		cfg,
		xyz[arg.src[0]],
		qzsign, abs(q[2]), qysign, abs(q[1]), qxsign, abs(q[0]),
		arg.src[1]);

	std::cout << __func__ << ": filen=" << filen << std::endl;

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
struct RhoPionFFSysParam
{
  RhoPionFFSysParam();
  RhoPionFFSysParam(XMLReader& xml_in, const std::string& path);

  struct Manage3PtParam
  {
    std::string  pattern;
    std::string  cache_file;
    std::string  cfg_file;
  } threept;

  struct Manage2PtParam
  {
    Array<std::string>   pattern_Z_pion;
    bool                 avg_2pt_func;
  } twopt;

  struct ManageEnergyParam
  {
    Real                 photon_Qsq;
    Array<std::string>   pattern_pion;
  } E;

  LatticeParam       lattice;

  string    Z_V_file_name;
  string    qsq_file_name;

  int       t_source_start;
  int       t_source_end;
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




// Rho-pi sys params
RhoPionFFSysParam::RhoPionFFSysParam(XMLReader& xml_in, const std::string& path) 
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

      read(top, "pattern_Z_pion", twopt.pattern_Z_pion);
      read(top, "avg_2pt_func", twopt.avg_2pt_func);
    }

    {
      XMLReader top(paramtop, "Energy");

      //  read(top, "photon_mass", E.photon_mass);
      read( top, "LSZ_photon_Qsq", E.photon_Qsq);

      read(top, "pattern_pion", E.pattern_pion);
    }

    read(paramtop, "LatticeParam", lattice);
    read(paramtop, "Z_V_file_name", Z_V_file_name);
    read(paramtop, "qsq_file_name", qsq_file_name);
    read(paramtop, "t_source_start", t_source_start);
    read(paramtop, "t_source_end", t_source_end);
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
class RhoPionFFSys: public LLSqComponent<EnsemReal>
{
public:
  //! destructor to help with cleanup;
  RhoPionFFSys(const RhoPionFFSysParam& param_);

  //! destructor to help with cleanup;
  ~RhoPionFFSys() {}

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

  //! Number of levels at sink
  int nLevelsF() const {return n_levels_f;}

  //! Set number of levels at source
  void setNLevelsF(int n);
  
  //! Return the dimensions of all the indices
  const ArrayInt& size() const {return sz;}

  //! Energy of initial state at momentum p
  EnsemReal getMassI() const;

  //! Energy of final state at momentum p
  EnsemReal& getMassF(int m) const;

protected:
  //! Convenience - make a threept arg
  ThreePtArg threePtArg(const PiPf& pi_pf, int g, int nu, int t_i) const;

  //! Convenience - make a twopt arg
  TwoPtArg twoPtArg(const ArrayInt& p) const;

  //! Energy of photon at momentum p
  EnsemReal energy_i(const ArrayInt& p) const;

  //! Energy of pi at momentum p
  EnsemReal& energy_f(int level, const ArrayInt& p) const;

  //! Snarf the threept function for a given source time
  EnsemReal getThreePt(int p, int mu, int nu, int t_i, int t) const;

  //! Compute matrix element piece
  EnsemReal matelem(int m, const ArrayInt& p_f, const ArrayInt& p_i, 
		    int mu, int nu) const;

  //! Compute propagation factors
  EnsemReal propfact(int m, const ArrayInt& p_f, const ArrayInt& p_i, int t) const;

protected:
  //! Hide default onstructor
  RhoPionFFSys(); 

private:
  RhoPionFFSysParam  param;

  mutable Manage3PtFuncBB      threept;
  mutable Array<ManageEnergy>  E_pion;
  mutable Array<ManageZAmp>    Z_pion;

  EnsemReal    Z_V;

  int              nUnknown;
  int              maxNUnknown;
  int              n_levels_f;
  int              nbin;
  int              t_fit_i;
  int              t_fit_f;
  Array<PiPf>      pi_pf;
  ArrayInt         sz;
};



// Constructor
RhoPionFFSys::RhoPionFFSys(const RhoPionFFSysParam& param_) : 
  param(param_), 
  threept(new RhoPion3PtFunc(param.threept.pattern), param.threept.cache_file, param.threept.cfg_file)
{
  t_fit_i = param.fit_range[0][0];   //! start of fit range
  t_fit_f = param.fit_range[0][1];   //! end of fit range

  // Read Z_V
  read(param.Z_V_file_name, Z_V);
  nbin = Z_V.size();
  cout << __func__ << ": found nbins=" << nbin << endl;

  // Read energies
  E_pion.resize(param.E.pattern_pion.size());
  for(int i=0; i < param.E.pattern_pion.size(); ++i)
    E_pion[i].create(new LocalEnergyFunc(param.E.pattern_pion[i], "pion"));

  // Read amplitudes
  Z_pion.resize(param.twopt.pattern_Z_pion.size());
  for(int i=0; i < param.twopt.pattern_Z_pion.size(); ++i)
    Z_pion[i].create(new LocalZFunc(param.twopt.pattern_Z_pion[i], "pion"), param.twopt.avg_2pt_func);

  // Check array sizes - number of exponentials
  if (Z_pion.size() !=  E_pion.size())
  {
    cerr << __func__ << ": inconsistent array size for E and Z" << endl;
    exit(1);
  }

  // Check ground state energies  
  ArrayInt zero(Nd-1);
  zero = 0;
  if (E_pion[0][zero].size() != nbins())
  {
    cerr << __func__ << ": energies inconsistent size with Z_V" << endl;
    exit(1);
  }

  n_levels_f = E_pion.size();
  maxNUnknown = 2 * E_pion.size();
  nUnknown    = 2;
  
  cout << __func__ << ": number of energy levels at sink   = " << nLevelsF() << endl;
  cout << __func__ << ": maximum number of unknowns = " << maxNUnknown << endl;
}


//! Convenience - make a threept arg
ThreePtArg RhoPionFFSys::threePtArg(const PiPf& pi_pf, int g, int nu, int t_i) const
{
  ThreePtArg  arg;
  arg.pi_pf  = pi_pf;
  arg.g      = g;
  arg.src.resize(2);
  arg.src[0] = nu;
  arg.src[1] = t_i;

//  std::cout << "3pt[]:  g= " << arg.g 
//	    << "   p_f= " << arg.pi_pf.p_f 
//	    << "   p_i= " << arg.pi_pf.p_i
//	    << "   nu= "  << nu
//	    << "   t_i= " << t_i
//	    << std::endl;

  return arg;
}

  //! Convenience - make a twopt arg
TwoPtArg RhoPionFFSys::twoPtArg(const ArrayInt& p) const
{
  TwoPtArg  arg;
  arg.p = p;
  return arg;
}



//! Set number of FF
void RhoPionFFSys::setNUnknowns(int n) 
{
  if (n > maxNUnknowns())
  {
    cerr << __func__ << ": too many unknowns\n";
    exit(1);
  }
  nUnknown = n;
}


// Reset Qsq to be used
void RhoPionFFSys::setQsq(const QsqVal& qsq)
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
  sz[0] = t_fit_f - t_fit_i + 1;   // fitting time range
  sz[1] = Nd;              // directions
  sz[2] = Nd-1;            // number of spatial directions
  sz[3] = pi_pf.size();    // number of momenta pairs
}


//! Reset time fit interval
void RhoPionFFSys::setTimeFit(int ti, int tf) 
{
  t_fit_i = ti; t_fit_f = tf;

  sz[0] = t_fit_f - t_fit_i + 1;   // reset time range
}


//! Get time fit interval
ArrayInt RhoPionFFSys::getTimeFit() const
{
  ArrayInt tt(2);
  tt[0] = t_fit_i;
  tt[1] = t_fit_f;

  return tt;
}

//! Get time fit interval
const Array<ArrayInt>& RhoPionFFSys::getFitRange() const
{
  return param.fit_range;
}


//! Set number of levels at sink
void RhoPionFFSys::setNLevelsF(int n)
{
  if (n > E_pion.size())
  {
    cerr << __func__ << ": number of levels too high\n";
    exit(1);
  }
  n_levels_f = n;
}
  

//! Compute photon energy using dispersion relation
EnsemReal RhoPionFFSys::energy_i(const ArrayInt& p) const
{
//  cout << __func__ << endl;

  //getMassI() in this case will return the Qsq value

  //NB The fitter doesn't return the right value for the other Qsq - need to modify this eventually

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
	 e =  photonenergyDisp(getMassI(), p, param.lattice);
        }
		   else
        {
	      e = Real(0.0);
        }


	 return e ;

  
}


//! Compute or return energy
EnsemReal& RhoPionFFSys::energy_f(int level, const ArrayInt& p0) const
{
//  cout << __func__ << endl;

  ArrayInt p = canonicalOrder(p0);
  ArrayInt zero(Nd-1);
  zero = 0;

#if 0
  if (p != zero && ! E_pion[level].exist(p))
    E_pion[level].insert(p, energyDisp(E_pion[level][zero], p, param.lattice));
#endif

  return E_pion[level][p];
}


//! Energy of initial state at momentum p
EnsemReal RhoPionFFSys::getMassI() const
{
  //actually the Qsq value for the LSZ photon
  EnsemReal Qsq;
  Qsq.resize(nbins());
  Qsq = param.E.photon_Qsq;
  return Qsq;
}


//! Energy of final state at momentum p
EnsemReal& RhoPionFFSys::getMassF(int m) const
{
  ArrayInt zero(Nd-1);
  zero = 0;
  return E_pion[m][zero];
}


//! Compute matrix element piece
EnsemReal RhoPionFFSys::matelem(int m, const ArrayInt& p_f, const ArrayInt& p_i, 
				int mu, int nu) const
{
  if (nu >= Nd-1)
  {
    cerr << __func__ << ": currently only support mu=0..3, nu<3\n";
    exit(1);
  }
  
  const double xi  = param.lattice.aniso.xi;
  //! IMPLEMENT SCALING OF MOMENTA BY c - NB will change the Qsq values too
  const double csq = param.lattice.aniso.c_sq;
  ArrayDouble pp_i = sqrt(csq) * contMom(p_i,param.lattice.latt_size[0]) / xi;
  ArrayDouble pp_f = sqrt(csq) * contMom(p_f,param.lattice.latt_size[0]) / xi;

  ArrayInt zero(Nd-1);
  zero = 0;

  const EnsemReal& m_f = getMassF(m);
  const EnsemReal  m_i = getMassI();
  const EnsemReal& E_f = energy_f(m,p_f);
  const EnsemReal& E_i = energy_i(p_i);
  EnsemReal        rat = Real( sqrt( xi*xi*xi) ) / (m_f);   //This is the normalisation I've been using to compare to data, includes the scaling for the anisotropy

  // Transition matrix element
  EnsemReal f;
  f.resize(nbins());
  f = 0;

  // Really dumb method
  TensorReal anti = antiSymTensor<4>();
  Array<int> n1(4); 

  for(int alpha=0; alpha < Nd; ++alpha)
    for(int delta=0; delta < Nd; ++delta)
    {
      n1[0] = mu; n1[1] = alpha; n1[2] = delta; n1[3] = nu;

      Real a = peekObs(anti,n1);

      // printf("%s: mu=%d alpha=%d delta=%d nu=%d:     a= %d\n",
      //   __func__,mu,alpha,delta,nu,int(toFloat(a)));

      if (toFloat(a) != 0)
      {
	if (alpha == Nd-1)
	{
	  f += a * rat*E_f*Real(pp_i[delta]);
	}
	else if (delta == Nd-1)
	{
	  f += a * rat*Real(pp_f[alpha])*E_i;
	}
	else 
	{
	  f += a * rat*Real(pp_f[alpha])*Real(pp_i[delta]);
	}

	//printf("%s: mu=%d alpha=%d delta=%d nu=%d: pf=%d,%d,%d, pi=%d,%d,%d: a=%d f=%f\n",
	// __func__,mu,alpha,delta,nu,
	//  p_f[0],p_f[1],p_f[2],p_i[0],p_i[1],p_i[2],
	// int(toFloat(a)),toFloat(mean(f)));

      }
    }

 

  return f;
}



//! Compute propagation factors
EnsemReal RhoPionFFSys::propfact(int m, const ArrayInt& p_f, const ArrayInt& p_i, int t) const
{
  int T_f = param.t_sink;

  EnsemReal g;
  g  = Z_pion[m][twoPtArg(p_f)] * exp(-Real(abs(T_f - t))*energy_f(m,p_f));
  g /= Real(2) * energy_f(m,p_f);

  return g;
}


// Given source coordinates, return the corresponding matrix element
EnsemReal RhoPionFFSys::mat(const ArrayInt& ind, int nFF) const
{
//  cout << __func__ << endl;

  //
  // nFF ==  0 ground state spatial V
  // nFF ==  1 ground state temporal V
  // nFF ==  2 first excited state spatial V
  // nFF ==  3 first excited state temporal V
  // nFF = <spatial/temporal V> + 2*m
  //
  int t  = ind[0] + t_fit_i;
  int mu = ind[1];
  int nu = ind[2];
  int p  = ind[3];

  ArrayInt p_i = pi_pf[p].p_i;
  ArrayInt p_f = pi_pf[p].p_f;

  int m = nFF >> 1; // Determine sink energy levels

  EnsemReal f;
  f.resize(nbins());
  f = 0;

  // NOTE: should eventually support cross term in transition FF with additional propfacts

  // spatial V
  if ((nFF & 1) == 0 && (mu < Nd-1))
  {
    f = matelem(m, p_f, p_i, mu, nu) * propfact(m, p_f, p_i, t);
  }

  // temporal V
  if ((nFF & 1) == 1 && (mu == Nd-1))
  {
    f = matelem(m, p_f, p_i, mu, nu) * propfact(m, p_f, p_i, t);
  }

  // Final kinetic factor
  f /= Z_V*Z_V;

//  cout << "exiting mat" << endl;

  return f;
}


//! Snarf the threept function for a given source time
EnsemReal RhoPionFFSys::getThreePt(int p, int mu, int nu, int t_i, int t) const
{
  //
  // Snarf the threept
  //
  EnsemVectorComplexF& three = threept[threePtArg(pi_pf[p], 1 << mu, nu, t_i)];
  EnsemReal thr;

  // Compensate for minus sign in rho_y from sequential source construction
  // This sign is independent of mu or momenta
  int sg = (nu == 1)? -1 : 1;

  // Extract appropriate part of 3pt. Fold in rho_y sign above.
  // signs are such that a minus sign in the three-point defn has cancelled with a minus sign to correct for the euclidean antisym four-tensor 
  if (mu == Nd-1)
    thr =  Real(sg)*real(peekObs(three,t));
  else
    thr =  Real(sg)*imag(peekObs(three,t));

  return thr;
}


// Given source coordinates, return the corresponding rhs element
EnsemReal RhoPionFFSys::rhs(const ArrayInt& ind) const
{
//  cout << __func__ << endl;

  int t  = ind[0] + t_fit_i;
  int mu = ind[1];
  int nu = ind[2];
  int p  = ind[3];

  if (nu < 0 || nu >= Nd-1)
  {
    cerr << __func__ << ": nu out of bounds\n";
    exit(1);
  }

  EnsemReal f;
  f.resize(nbins());
  f = 0;

  ArrayInt p_i = pi_pf[p].p_i;
  ArrayInt p_f = pi_pf[p].p_f;

  // Integrate over source time
  for(int t_i=param.t_source_start; t_i <= param.t_source_end; ++t_i)
  {
    f += getThreePt(p, mu, nu, t_i, t) * exp(-Real(t_i - t)*energy_i(p_i));
  }

//  cout << "exiting mat" << endl;

  return f;
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
  RhoPionFFSysParam linParam(xml_in, "/Pion2Photon");

  XMLFileWriter xml_out(argv[2]);
  push(xml_out, "Pion2Photon");

  cout << "construct linear system" << endl;
  RhoPionFFSys sys(linParam);

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

    //this is not the mass anymore - its the Qsq of the LSZ photon!
    qsqParam.mass_i         = toDouble(mean(sys.getMassI()));
  qsqParam.mass_f         = toDouble(mean(sys.getMassF(0)));
  
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
    const double mass_f = toDouble(mean(sys.getMassF(0)));
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

  pop(xml_out);  // RhoPionFF
  xml_out.close();
}
