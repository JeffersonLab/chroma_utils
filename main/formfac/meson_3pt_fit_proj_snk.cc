// $Id: meson_3pt_fit_proj_snk.cc,v 2.3 2009/03/25 02:24:15 edwards Exp $
/*! \file
 * \brief Meson-Meson form-factor code using fitting method
 */

#include "adat/handle.h"
#include "adat/map_obj.h"
#include "io/adat_xmlio.h"
#include "io/adat_io.h"
#include "io/adat_xml_group_reader.h"
#include "ensem/ensem.h"
#include "formfac/formfac_manage_2pt.h"
#include "formfac/formfac_manage_E.h"
#include "formfac/formfac_solver.h"
#include "formfac/formfac_solver_driver.h"
#include "formfac/formfac_manage_factory.h"
#include "formfac/formfac_manage_aggregate.h"
#include "current_ops_factory.h"
#include "operators_factory.h"
#include "matrix_elems_factory.h"
#include <list>
#include <iostream>

using namespace ADAT;
using namespace ENSEM;
using namespace FF;
using namespace std;

namespace FF
{
  //! Print an array
  std::ostream& operator<<(std::ostream& s, const Array<int>& d)
  {
    s << d[0];
    for(int i=1; i < d.size(); ++i)
      s << " " << d[i];

    return s;
  }


  //! Print an array
  std::ostream& operator<<(std::ostream& s, const Array<double>& d)
  {
    s << d[0];
    for(int i=1; i < d.size(); ++i)
      s << " " << d[i];

    return s;
  }


  //-----------------------------------------------------------------------------------
  //! Params for linear system
  struct MesonFFSysParam
  {
    MesonFFSysParam();
    MesonFFSysParam(XMLReader& xml_in, const std::string& path);

    /*! 
     * Being bold here - instead of an array of ground+excited states, just 
     * list them all as separate states. This allows only excited state fits with
     * no ground state!
     */
    struct SourceSinkOperators_t
    {
      //! All these operators share the same energy level
      struct Operator_t
      {
	GroupXML_t       op;            /*!< xml holding operator group */
	string           smear;         /*!< Smearing label to be used */
      };

      GroupXML_t         energy;        /*!< xml holding E group */
      GroupXML_t         amp;           /*!< xml holding Z group */
      string             quantum_id;    /*!< quantum number label for overlap-state */
      Array<Operator_t>  operators;     /*!< Operators used in determining this state */
    };

    SourceSinkOperators_t  source_operators;
    SourceSinkOperators_t  sink_operators;

    struct InsertOp_t
    {
      int          src_level;     /*!< The source should be at this energy level */
      int          snk_level;     /*!< The sink should be at this energy level */
      GroupXML_t   matrix_elem;   /*!< xml holding matrix element */
      GroupXML_t   current_op;    /*!< xml holding matrix element */
    };

    Array<InsertOp_t> insertion_operators;

    struct AdditionalInfo_t
    {
      int          dt;           /*!< Source-sink separation */
      int          quark;        /*!< Quark line number for the insertion */
      string       mass;         /*!< Mass label used for the 3-pt correlator */
      string       ensemble;     /*!< Ensemble label used in the dbs */
    };

    AdditionalInfo_t   db_info;
    LatticeParam       lattice;

    string             qsq_file_name;

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

    std::vector<TimeFitRange_t>  fit_ranges;
  };



  // Source/sink operators
  void read(XMLReader& xml_in, const std::string& path, MesonFFSysParam::SourceSinkOperators_t::Operator_t& param) 
  {
    XMLReader paramtop(xml_in, path);

    read(paramtop, "smear_label", param.smear);
    param.op  = ADATXML::readXMLGroup(paramtop, "Operator", "OperatorType");
  }


  // Source/sink operators
  void read(XMLReader& xml_in, const std::string& path, MesonFFSysParam::SourceSinkOperators_t& param) 
  {
    XMLReader paramtop(xml_in, path);

    read(paramtop, "Operators", param.operators);
    read(paramtop, "quantum_id", param.quantum_id);
    param.energy = ADATXML::readXMLGroup(paramtop, "Energy", "EnergyType");
    param.amp = ADATXML::readXMLGroup(paramtop, "Amplitude", "AmpType");
  }


  // Source/sink operators
  void read(XMLReader& xml_in, const std::string& path, MesonFFSysParam::InsertOp_t& param) 
  {
    XMLReader paramtop(xml_in, path);

    read(paramtop, "source_level", param.src_level);
    read(paramtop, "sink_level", param.snk_level);
    param.matrix_elem = ADATXML::readXMLGroup(paramtop, "MatrixElement", "MatrixElemType");
    param.current_op  = ADATXML::readXMLGroup(paramtop, "CurrentOperator", "CurrentOperatorType");
  }


  // Additional info used for the dbs
  void read(XMLReader& xml_in, const std::string& path, MesonFFSysParam::AdditionalInfo_t& param) 
  {
    XMLReader paramtop(xml_in, path);

    read(paramtop, "dt", param.dt);
    read(paramtop, "quark", param.quark);
    read(paramtop, "mass", param.mass);
    read(paramtop, "ensemble", param.ensemble);
  }



  // Meson-meson sys params
  MesonFFSysParam::MesonFFSysParam(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);

    read(paramtop, "SourceOperators", source_operators);
    read(paramtop, "SinkOperators", sink_operators);
    read(paramtop, "InsertionOperators", insertion_operators);
    read(paramtop, "AdditionalInfo", db_info);
    read(paramtop, "LatticeParam", lattice);
    read(paramtop, "qsq_file_name", qsq_file_name);
    read(paramtop, "fit_ranges", fit_ranges);

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


  //-----------------------------------------------------------------------------------
  //! Hold operator and wavefunctions
  struct OperWaveFuncKey_t
  {
    string    name;    /*!< Operator name */
    string    smear;   /*!< Operator smearing label */
  };

  //! Hold operator and wavefunctions
  struct OperWaveFuncVal_t
  {
    Handle< Operator<Complex> >       op;
    Handle< WaveFunction<Complex> >   wvf;
  };



  //----------------------------------------------------------------------------
  // Support for hashes
  void read(BinaryReader& bin, OperWaveFuncKey_t& param)
  {
    readDesc(bin, param.name);
    readDesc(bin, param.smear);
  }

  void write(BinaryWriter& bin, const OperWaveFuncKey_t& param)
  {
    writeDesc(bin, param.name);
    writeDesc(bin, param.smear);
  }


  //----------------------------------------------------------------------------------
  // Concatenate two vectors. 
  // Am I missing something? I thought this should be available already.
  template<typename T>
  inline std::vector<T> concat(const std::vector<T>& s1, const std::vector<T>& s2)
  {
    std::vector<T> dest;

    // Loop over each source term
    for(int i=0; i < s1.size(); ++i)
    {
      dest.push_back(s1[i]);
    }

    // Loop over each source term
    for(int i=0; i < s2.size(); ++i)
    {
      dest.push_back(s2[i]);
    }
    return dest;
  }

  //-----------------------------------------------------------------------------------
  //! System for linear least squares
  /*! System here is rho-pi FF */
  class MesonFFSys: public LLSqComponent<EnsemReal>
  {
  public:
    //! Full constructor
    MesonFFSys(const MesonFFSysParam& param_);

    //! Destructor
    ~MesonFFSys() {}

    //! The names of the unknowns. These are the possible column indices
    std::vector<std::string> unknowns() const {return FF_names;}

    //! Set qsq to use
    void setQsq(const QsqVal& qsq);

    //! Reset time fit interval
    void setTimeFit(const TimeFitRange_t& range);

    //! Get time fit interval
    TimeFitRange_t getTimeFit() const;

    //! Get time fit interval
    std::vector<TimeFitRange_t> getFitRanges() const;

    //! Time extent
    int timeLen() const {return param.lattice.latt_size[param.lattice.decay_dir];}

    //! Return the list of the row indices
    std::vector<LLSqRow_t> rowList() const;

    //! Given row and column indices, return the corresponding matrix element
    EnsemReal mat(const LLSqRow_t& row, const std::string& FF_name) const;

    //! Given row index, return the corresponding rhs element
    EnsemReal rhs(const LLSqRow_t& row) const;

    //! Number of bins
    int nbins() const {return nbin;}

    //! Possibly save the solutions
//    void save(const Array<EnsemReal>& F, 
//	      const Real& qsq, int t_i, int t_f) const;

    //! Energy of initial state at zero momentum
    EnsemReal getMassI(int n) const;

    //! Energy of final state at zero momentum
    EnsemReal getMassF(int m) const;

  protected:
    //! Convenience - make a energy arg
    EnergyArg energyArg(int level, const ArrayInt& p) const;

    //! Convenience - make a amp arg
    AmpArg ampArg(int level, const ArrayInt& p, const string& name, const string& smear) const;

    //! Energy of initial state at momentum p
    EnsemReal energy_i(int level, const ArrayInt& p) const;

    //! Energy of final state at momentum p
    EnsemReal energy_f(int level, const ArrayInt& p) const;

    //! Amp of initial state at momentum p
    EnsemReal amp_i(int level, const ArrayInt& p, const string& name, const string& smear) const;

    //! Amp of final state at momentum p
    EnsemReal amp_f(int level, const ArrayInt& p, const string& name, const string& smear) const;

    //! Given source coordinates, return the corresponding matrix element
    Real matelem_wvf(const MesonMatrixElement& matrix_elem, 
		     const LLSqRow_t& row, 
		     const std::string& FF_name,
		     const Real& m_f, const Real& m_i) const;

    //! Given row and column indices, return the corresponding matrix element
    EnsemReal matelem(const MesonMatrixElement& matrix_elem, 
		      const LLSqRow_t& row, 
		      const std::string& FF_name,
		      int m, int n) const;

    //! Compute propagation factors
    EnsemReal propfact(const LLSqRow_t& row, int m, int n, int t) const;

    //! Initialize operators and wavefuncs
    void initSourceOperators();

    //! Initialize operators and wavefuncs
    void initSinkOperators();

    //! Initialize energy and amp map
    void initSourceEnergyAmpMap();

    //! Initialize energy and amp map
    void initSinkEnergyAmpMap();

    //! Initialize insertion
    void initInsertion();

    //! Determine the form-factor names
    void initFFNames();

    //! Some basic sanity checks
    void sanityChecks();

  protected:
    //! Hide default onstructor
    MesonFFSys(); 

  private:
    MesonFFSysParam  param;


    MapObject<OperWaveFuncKey_t,OperWaveFuncVal_t>    source_operators;
    MapObject<OperWaveFuncKey_t,OperWaveFuncVal_t>    sink_operators;

    Handle< ManageAmpFuncMap >          Z_src;
    Handle< ManageEnergyFuncMap >       E_src;
    Handle< ManageAmpFuncMap >          Z_snk;
    Handle< ManageEnergyFuncMap >       E_snk;
  

    struct InsertionOperator_t
    {
      int                           src_level;
      int                           snk_level;
      Handle<MesonMatrixElement>    matrix_elem;
      Handle<CurrentOperator>       current_op;
    };

    std::vector<InsertionOperator_t>    insertion_operators;
    std::vector<std::string>            FF_names;


    int                       nbin;
    TimeFitRange_t            fit_range;
    std::vector<PiPf>         pi_pf;
  };



  //-----------------------------------------------------------------------------------
  // Setup the insertion and current operators
  void MesonFFSys::initInsertion()
  {
    //
    // Construct insertion operator object
    //
    for(int i=0; i < param.insertion_operators.size(); ++i)
    {
      const MesonFFSysParam::InsertOp_t& op_xml = param.insertion_operators[i];
      InsertionOperator_t op;

      op.src_level = op_xml.src_level;
      op.snk_level = op_xml.snk_level;

      // Sanity check
//      if (op.src_level+1 > param.source_levels.size())
//      {
//	cerr << __func__ << ": source level for insertion greater than number of source energy levels\n";
//	exit(1);
//      }
//
//      if (op.snk_level+1 > param.sink_levels.size())
//      {
//	cerr << __func__ << ": sink level for insertion greater than number of sink energy levels\n";
//	exit(1);
//      }

      // Construct the matrix element
      try
      {
	std::istringstream  xml_s(op_xml.matrix_elem.xml);
	XMLReader  optop(xml_s);
	
	op.matrix_elem = TheMesonMatrixElementFactory::Instance().createObject(op_xml.matrix_elem.id,
									       optop,
									       op_xml.matrix_elem.path);
      }
      catch(const std::string& e) 
      {
	cerr << "Caught Exception reading matrix element XML: " << e << endl;
	exit(1);
      }
 
      // Construct the current operator
      try
      {
	std::istringstream  xml_s(op_xml.current_op.xml);
	XMLReader  optop(xml_s);
	
	op.current_op = TheCurrentOperatorFactory::Instance().createObject(op_xml.current_op.id,
									   optop,
									   op_xml.current_op.path);
      }
      catch(const std::string& e) 
      {
	cerr << "Caught Exception reading current operator XML: " << e << endl;
	exit(1);
      }

      // Now a valid entry, push it on the back of all the insertion operators
      insertion_operators.push_back(op);
    }
  }


  // Setup all the operators and wavefunctions
  void MesonFFSys::initSourceOperators()
  {
    //
    // Create the source operator and wave function overlaps
    //
    for(int src_op=0; src_op < param.source_operators.operators.size(); ++src_op)
    {
      const MesonFFSysParam::SourceSinkOperators_t::Operator_t& in_op = param.source_operators.operators[src_op];
      OperWaveFuncKey_t key;
      OperWaveFuncVal_t val;

      try
      {
	std::istringstream  xml_s(in_op.op.xml);
	XMLReader  optop(xml_s);
	
	val.op = TheMesonOperatorFactory::Instance().createObject(in_op.op.id,
								  optop,
								  in_op.op.path);
      }
      catch(const std::string& e) 
      {
	cerr << "Caught Exception source operator XML: " << e << endl;
	exit(1);
      }

      // Construct the wave functions, and find the one with the matching quantum id
      std::vector< Handle< WaveFunction<Complex> > > wvfs(val.op->overlaps());

      int fnd = 0;
      for(int i=0; i < wvfs.size(); ++i)
      {
	if (wvfs[i]->particleName() == param.source_operators.quantum_id)
	{
	  val.wvf = wvfs[i];
	  ++fnd;
	}
      }
      if (fnd != 1)
      {
	cerr << __func__ << ": did not find quantum_id= " << param.source_operators.quantum_id 
	     << "  in wavefunction list for source operator= " << in_op.op.id
	     << endl;
	exit(1);
      }

      // New entry
      key.name  = in_op.op.id;
      key.smear = in_op.smear;

      source_operators.insert(key,val);
    } // for src_op
  }


  // Setup all the operators and wavefunctions
  void MesonFFSys::initSinkOperators()
  {
    // At the moment, will only support one sink operator. 
    // The idea is that the sink will be projected, so there is only one state.
    // This probably needs to be generalized if sink is at non-zero momentum,
    // in which case the there might be additional overlaps for some p_f and
    // sink direction.
    if (param.sink_operators.operators.size() != 1)
    {
      cerr << __func__ << ": only support 1 sink state\n";
      exit(1);
    }

    //
    // Create the sink operator and wave function overlaps
    //
    for(int snk_op=0; snk_op < param.sink_operators.operators.size(); ++snk_op)
    {
      const MesonFFSysParam::SourceSinkOperators_t::Operator_t& in_op = param.sink_operators.operators[snk_op];
      OperWaveFuncKey_t key;
      OperWaveFuncVal_t val;

      try
      {
	std::istringstream  xml_s(in_op.op.xml);
	XMLReader  optop(xml_s);
	
	val.op = TheMesonOperatorFactory::Instance().createObject(in_op.op.id,
								  optop,
								  in_op.op.path);
      }
      catch(const std::string& e) 
      {
	cerr << "Caught Exception sink operator XML: " << e << endl;
	exit(1);
      }

      // Construct the wave functions, and find the one with the matching quantum id
      std::vector< Handle< WaveFunction<Complex> > > wvfs(val.op->overlaps());

      int fnd = 0;
      for(int i=0; i < wvfs.size(); ++i)
      {
	if (wvfs[i]->particleName() == param.sink_operators.quantum_id)
	{
	  val.wvf = wvfs[i];
	  ++fnd;
	}
      }
      if (fnd != 1)
      {
	cerr << __func__ << ": did not find quantum_id= " << param.sink_operators.quantum_id 
	     << "  in wavefunction list for sink operator= " << in_op.op.id
	     << endl;
	exit(1);
      }

      // New entry
      key.name  = in_op.op.id;
      key.smear = in_op.smear;

      sink_operators.insert(key,val);
    } // for snk_op
  }


  // Setup the energy and amplitude maps
  void MesonFFSys::initSourceEnergyAmpMap()
  {
    // Create the energy map
    try
    {
      std::istringstream  xml_s(param.source_operators.energy.xml);
      XMLReader  optop(xml_s);
	
      E_src = TheManageEnergyFuncMapFactory::Instance().createObject(param.source_operators.energy.id,
								     optop,
								     param.source_operators.energy.path);
    }
    catch(const std::string& e) 
    {
      cerr << "Caught Exception reading source energy XML: " << e << endl;
      exit(1);
    }

    // Create the amp map
    try
    {
      std::istringstream  xml_s(param.source_operators.amp.xml);
      XMLReader  optop(xml_s);
	
      Z_src = TheManageAmpFuncMapFactory::Instance().createObject(param.source_operators.amp.id,
								  optop,
								  param.source_operators.amp.path);
    }
    catch(const std::string& e) 
    {
      cerr << "Caught Exception reading source amp XML: " << e << endl;
      exit(1);
    }
  }


  // Setup the energy and amplitude maps
  void MesonFFSys::initSinkEnergyAmpMap()
  {
    // Create the energy map
    try
    {
      std::istringstream  xml_s(param.sink_operators.energy.xml);
      XMLReader  optop(xml_s);
	
      E_snk = TheManageEnergyFuncMapFactory::Instance().createObject(param.sink_operators.energy.id,
								     optop,
								     param.sink_operators.energy.path);
    }
    catch(const std::string& e) 
    {
      cerr << "Caught Exception reading sink energy XML: " << e << endl;
      exit(1);
    }

    // Create the amp map
    try
    {
      std::istringstream  xml_s(param.sink_operators.amp.xml);
      XMLReader  optop(xml_s);
	
      Z_snk = TheManageAmpFuncMapFactory::Instance().createObject(param.sink_operators.amp.id,
								  optop,
								  param.sink_operators.amp.path);
    }
    catch(const std::string& e) 
    {
      cerr << "Caught Exception reading sink amp XML: " << e << endl;
      exit(1);
    }
  }


  // Determine the form-factor names
  void MesonFFSys::initFFNames()
  {
    FF_names.clear();

    for(int i=0; i < insertion_operators.size(); ++i)
    {
      const MesonMatrixElement& matrix_elem = *(insertion_operators[i].matrix_elem);
      FF_names = concat(FF_names, matrix_elem.getFFNames());
    }
  }


  // Constructor
  MesonFFSys::MesonFFSys(const MesonFFSysParam& param_) : 
    param(param_)
  {
    initInsertion();
    initSourceOperators();
    initSinkOperators();
    initSourceEnergyAmpMap();
    initSinkEnergyAmpMap();
    initFFNames();

    sanityChecks();

    nbin = E_src->size();  // Should do more checks
  }


  // Sanity checks
  void MesonFFSys::sanityChecks()
  {
    // Check stuff like number of bins for each manager, etc.

  }

  //! Convenience - make a energy arg
  EnergyArg MesonFFSys::energyArg(int level, const ArrayInt& p) const
  {
    EnergyArg arg;
    arg.level = level;
    arg.mom = p;
    return arg;
  }


  //! Convenience - make an amp arg
  AmpArg MesonFFSys::ampArg(int level, const ArrayInt& p, const string& name, const string& smear) const
  {
    AmpArg arg;
    arg.level = level;
    arg.mom = p;
    arg.name = name;
    arg.smear = smear;
    return arg;
  }



  // Reset Qsq to be used
  void MesonFFSys::setQsq(const QsqVal& qsq)
  {
    int n = 0;
    for(list<PiPf>::const_iterator nqq=qsq.pi_pf.begin(); 
	nqq != qsq.pi_pf.end(); 
	++nqq,++n)
    {
      pi_pf.push_back(*nqq);
    }
  }


  //! Reset time fit interval
  void MesonFFSys::setTimeFit(const TimeFitRange_t& range_) 
  {
    fit_range = range_;
  }


  //! Get time fit interval
  TimeFitRange_t MesonFFSys::getTimeFit() const
  {
    return fit_range;
  }


  //! Get time fit interval
  std::vector<TimeFitRange_t> MesonFFSys::getFitRanges() const
  {
    return param.fit_ranges;
  }


  //! Return the list of the row indices
  std::vector<LLSqRow_t> MesonFFSys::rowList() const
  {
    std::vector<LLSqRow_t> out;

    // Not sure I want to do this
    std::vector<OperWaveFuncKey_t> src_keys(source_operators.keys());
    std::vector<OperWaveFuncKey_t> snk_keys(sink_operators.keys());

    for(int t=fit_range.t_i; t <= fit_range.t_f; ++t)
    {
      for(int p=0; p < pi_pf.size(); ++p)
      {
	for(int ins=0; ins < insertion_operators.size(); ++ins)
	{
	  
	  for(int srck=0; srck < src_keys.size(); ++srck)
	  {
	    for(int snkk=0; snkk < snk_keys.size(); ++snkk)
	    {
	      // Build the key
	      LLSqRow_t row;

	      row.t           = t;
	      row.pi_pf       = pi_pf[p];
	      row.insert_lorentz = 0;   // FIXME

	      row.dt          = param.db_info.dt;
	      row.quark       = param.db_info.quark;

	      row.src_name    = src_keys[srck].name;
	      row.src_smear   = src_keys[srck].smear;
	      row.src_lorentz = 0;    // FIXME
	      row.src_spin    = -1;      // This comes from the 3pt meson code

	      row.snk_name    = snk_keys[snkk].name;
	      row.snk_smear   = snk_keys[snkk].smear;
	      row.snk_lorentz = 0;    // FIXME
	      row.snk_spin    = -1;      // This comes from the 3pt meson code

	      row.mass        = param.db_info.mass;
	      row.ensemble    = param.db_info.ensemble;

	      // Add this new entry to the list
	      out.push_back(row);
	    } // snkk
	  } // srck
	} // ins
      } // p 
    } // t

    
    return out;
  }

  //! Compute or return energy
  EnsemReal MesonFFSys::energy_i(int level, const ArrayInt& p0) const
  {
//  cout << __func__ << endl;

    return (*E_src)[energyArg(level, canonicalOrder(p0))];
  }


  //! Compute or return energy
  EnsemReal MesonFFSys::energy_f(int level, const ArrayInt& p0) const
  {
//  cout << __func__ << endl;

    return (*E_snk)[energyArg(level, canonicalOrder(p0))];
  }


  //! Compute or return energy
  EnsemReal MesonFFSys::amp_i(int level, const ArrayInt& p0, const string& name, const string& smear) const
  {
//  cout << __func__ << endl;

    return (*Z_src)[ampArg(level, canonicalOrder(p0), name, smear)];
  }


  //! Compute or return energy
  EnsemReal MesonFFSys::amp_f(int level, const ArrayInt& p0, const string& name, const string& smear) const
  {
//  cout << __func__ << endl;

    return (*Z_snk)[ampArg(level, canonicalOrder(p0), name, smear)];
  }


  //! Energy of initial state at momentum p
  EnsemReal MesonFFSys::getMassI(int level) const
  {
    ArrayInt zero(Nd-1);
    zero = 0;
    return energy_i(level, zero);
  }


  //! Energy of final state at momentum p
  EnsemReal MesonFFSys::getMassF(int level) const
  {
    ArrayInt zero(Nd-1);
    zero = 0;
    return energy_f(level, zero);
  }

 
  //! Compute propagation factors
  EnsemReal MesonFFSys::propfact(const LLSqRow_t& row, int m, int n, int t) const
  {
    int T_i = 0;
    int T_f = row.dt;
    int T   = param.lattice.latt_size[param.lattice.decay_dir];

    const ArrayInt& p_i = row.pi_pf.p_i;
    const ArrayInt& p_f = row.pi_pf.p_f;

    EnsemReal g, g_i, g_f;
    // now we've folded over the second half of the lattice we simplify the images (in fact they should be negligable)
    g_i = exp(-Real(t - T_i)*energy_i(n,p_i));       
    g_f = exp(-Real(abs(T_f - t))*energy_f(m,p_f));  
    g   = amp_f(m,p_f,row.snk_name,row.snk_smear) * amp_i(n,p_i,row.src_name,row.src_smear) * g_i * g_f;
    g  /= Real(4) * energy_f(m,p_f) * energy_i(n,p_i);

    return g;
  }


  // Given source coordinates, return the corresponding matrix element
  Real MesonFFSys::matelem_wvf(const MesonMatrixElement& matrix_elem, 
			       const LLSqRow_t& row, 
			       const std::string& FF_name,
			       const Real& m_f, const Real& m_i) const
  {
    cout << __func__ << ": entering\n";
    
    OperWaveFuncKey_t src_key;
    src_key.name  = row.src_name;
    src_key.smear = row.src_smear;
    const WaveFunction<Complex>& src_wave = *(source_operators[src_key].wvf);

    OperWaveFuncKey_t snk_key;
    snk_key.name  = row.snk_name;
    snk_key.smear = row.snk_smear;
    const WaveFunction<Complex>& snk_wave = *(sink_operators[snk_key].wvf);

    const ArrayInt& p_i = row.pi_pf.p_i;
    const ArrayInt& p_f = row.pi_pf.p_f;

    if (src_wave.numPolar() != matrix_elem.numSourcePolar())
    {
      std::cerr << __func__ << ": source polarizations do not match\n";
      exit(1);
    }

    int src_dir;
    if (row.src_lorentz.size() == 0)
      src_dir = -1;
    else if (row.src_lorentz.size() == 1)
      src_dir = row.src_lorentz[0];
    else
    {
      std::cerr << __func__ << ": src_lorentz size not 0 or 1\n";
      exit(1);
    }

    if (snk_wave.numPolar() != matrix_elem.numSinkPolar())
    {
      std::cerr << __func__ << ": source polarizations do not match\n";
      exit(1);
    }

    int snk_dir;
    if (row.snk_lorentz.size() == 0)
      snk_dir = -1;
    else if (row.snk_lorentz.size() == 1)
      snk_dir = row.snk_lorentz[0];
    else
    {
      std::cerr << __func__ << ": snk_lorentz size not 0 or 1\n";
      exit(1);
    }

    Array<Complex> wvf_snk(snk_wave.numPolar());
    Array<Complex> wvf_src(src_wave.numPolar());

    for(int r_f = 0; r_f < snk_wave.numPolar(); ++r_f)
    {
      Array<Complex> wvf = snk_wave(m_f, p_f, snk_dir, r_f);

      if (wvf.size() != 1)
      {
	std::cerr << __func__ << ": sink wavefunction size not equal to 1, cannot deal with it\n";
	exit(1);
      }

      wvf_snk[r_f] = wvf[0];
    }

    for(int r_i = 0; r_i < src_wave.numPolar(); ++r_i)
    {
      Array<Complex> wvf = src_wave(m_i, p_i, src_dir, r_i);

      if (wvf.size() != 1)
      {
	std::cerr << __func__ << ": source wavefunction size not equal to 1, cannot deal with it\n";
	exit(1);
      }

      wvf_src[r_i] = wvf[0];
    }

    Real f = zero;
    for(int r_f = 0; r_f < snk_wave.numPolar(); ++r_f)
    {
      for(int r_i = 0; r_i < src_wave.numPolar(); ++r_i)
      {
	MatElemRes_t<Complex,Complex> res = matrix_elem(FF_name, m_f, p_f, m_i, p_i, 
							row.insert_lorentz, r_f, r_i);

	Complex cf = wvf_snk[r_f] * res.result * wvf_src[r_i];
	f += real(cf);
      }
    }

    cout << __func__ << ": exiting\n";

    return f;
  }


  // Given source coordinates, return the corresponding matrix element
  EnsemReal MesonFFSys::matelem(const MesonMatrixElement& matrix_elem, 
				const LLSqRow_t& row, 
				const std::string& FF_name,
				int m, int n) const
  {
    cout << __func__ << ": entering\n";

    EnsemReal mass_f = rescaleEnsemDown(getMassF(m));
    EnsemReal mass_i = rescaleEnsemDown(getMassI(n));

    EnsemReal f;
    f.resize(nbins());

    for(int bin=0; bin < nbins(); ++bin)
    {
      Real m_f = peekEnsem(mass_f, bin);
      Real m_i = peekEnsem(mass_i, bin);

      // NEED TO CHECK REALITY OF MAT_ELEM
      Real ff = matelem_wvf(matrix_elem, row, FF_name, m_f, m_i);
      
      pokeEnsem(f, ff, bin);
    }

    cout << __func__ << ": exiting\n";

    return rescaleEnsemUp(f);
  }


  // Given source coordinates, return the corresponding matrix element
  EnsemReal MesonFFSys::mat(const LLSqRow_t& row, const std::string& FF_name) const
  {
    cout << __func__ << ": entering\n";

    const MesonMatrixElement& matrix_elem = *(insertion_operators[0].matrix_elem);

    int t = row.t;
    const ArrayInt& p_i = row.pi_pf.p_i;
    const ArrayInt& p_f = row.pi_pf.p_f;
    int m = 0;
    int n = 0;

    EnsemReal f = matelem(matrix_elem, row, FF_name, m, n) * propfact(row, m, n, t);

    cout << __func__ << ": exiting\n";

    return f;
  }


  // Given source coordinates, return the corresponding rhs element
  EnsemReal MesonFFSys::rhs(const LLSqRow_t& row) const
  {
    cout << __func__ << ": entering\n";

    const CurrentOperator& current_op = *(insertion_operators[0].current_op);

    EnsemReal thr = current_op(row);

    cout << __func__ << ": exiting\n";

    return thr;
  }


#if 0
  // Setup all the operators and wavefunctions
  void MesonFFSys::setup()
  {
    //
    // Create the source operator and wave function overlaps
    //
    source_operators.resize(param.source_operators.size());
    for(int src_op=0; src_op < param.source_operators.size(); ++src_op)
    {
      string source_id = param.source_operators[src_op].operator_id;
      cout << "Source operator = " << source_id << endl;

//    const string operator_path = "/FileTest";
      source_operators[src_op].op = TheMesonOperatorFactory::Instance().createObject(source_id,
										     xml_in,
										     operator_path);

      // Construct the wave functions
      source_operators[src_op].wvfs = source_operators[src_op].op->overlaps();
    }

    //
    // Create the sink operator and wave function overlaps
    //
    sink_operators.resize(param.sink_operators.size());
    for(int snk_op=0; snk_op < param.sink_operators.size(); ++snk_op)
    {
      string sink_id = param.sink_operators[snk_op].operator_id;
      cout << "Sink operator = " << sink_id << endl;

//    const string operator_path = "/FileTest";
      sink_operators[snk_op].op = TheMesonOperatorFactory::Instance().createObject(sink_id,
										   xml_in,
										   operator_path);

      // Construct the wave functions
      sink_operators[src_op].wvfs = sink_operators[src_op].op->overlaps();
    }


    //
    // Loop over source and sink wavefunctions, find the name (quantum number)
    // of the state and find allowed transition matrix elements
    //
    for(int src_op=0; src_op < param.source_operators.size(); ++src_op)
    {
      string source_op_name = source_operators[src_op].op->operatorName();
      cout << "source_op_name = " << source_op_name << endl;

      for(int snk_op=0; snk_op < param.sink_operators.size(); ++snk_op)
      {
	string sink_op_name = sink_operators[src_op].op->operatorName();
	cout << "sink_op_name = " << sink_op_name << endl;

	typedef list< Handle< WaveFunction<Complex> > > ListType_t;

	for(ListType_t::const_iterator src_wvf=param.source_operators[src_op].wvfs.begin();
	    src_wvf != param.source_operators[src_op].wvfs.end(); 
	    ++src_wvf)
	{
	  const WaveFunction<Complex>& source_wvf = *(*src_wvf);
	  string source_particle = source_wvf.particleName();
	  cout << "source_particle = " << source_particle << endl;

	  for(ListType_t::const_iterator src_wvf=param.sink_operators[snk_op].wvfs.begin();
	      snk_wvf != param.sink_operators[snk_op].wvfs.end(); 
	      ++snk_wvf)
	  {
	    const WaveFunction<Complex>& sink_wvf = *(*snk_wvf);
	    string sink_particle = sink_wvf.particleName();
	    cout << "sink_particle = " << sink_particle << endl;
    
	    // Use if the quantum numbers match matrix element
	    string mat_elem_id = matrixElementId(source_particle, sink_particle);
	    if (! useMatrixElementP(source_particle, sink_particle))
	    {
	      cout << __func__ << ": no matrix element for key = " << mat_elem_id << endl;
	      continue;
	    }

	    // We know the matrix element exist and intend to use it
	    Handle<MesonMatrixElement> mes_mat_elem(
	      TheMesonOperatorFactory::Instance().createObject(mat_elem_id,
							       xml_in,
							       operator_path));

	  }
	}
      }
    }
  }
#endif


#if 0
  //! Possibly save the solutions
  void MesonFFSys::save(const Array<EnsemReal>& F, 
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
      //   mkdir(dir, 0755);
      string makedir = string("mkdir -p ") + string(dir);
      system( makedir.c_str() );

      for(int n=0; n < nUnknowns(); ++n)
      {
	char ss[1024];
	EnsemReal f = F[n];

	sprintf(ss, param.saveSoln.pattern.c_str(),
		n, t_i, t_f);

	string filename = string(dir) + '/' + string(ss);
      
	write(filename, f);
      }
    }
  }
#endif


} // namespace FF


//using namespace ADATXML;
using namespace ENSEM;
using namespace FF;
using namespace std;


// Linkage hack
bool linkage()
{
  bool success = true;

  // Register all factories
  success &= FormfacManageEnergyFuncEnv::registerAll();
  success &= FormfacManageAmpFuncEnv::registerAll();
  success &= CurrentOperatorEnv::registerAll();
  success &= MesonOperatorEnv::registerAll();
  success &= MesonMatrixElementEnv::registerAll();

  return success;
}


int main(int argc, char *argv[])
{
  // Register all factories
  linkage();
  
  if (argc != 3)
  {
    cerr << "Usage:"
	 << argv[0] << ":  <input xml file>   <output xml file>" << endl;
    exit(1);
  }

  XMLReader xml_in(argv[1]);
  MesonFFSysParam linParam(xml_in, "/MesonFFFit");

  XMLFileWriter xml_out(argv[2]);
  push(xml_out, "MesonFFFit");

  try
  {
    cout << "construct linear system" << endl;
    MesonFFSys sys(linParam);

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

    double cov_tol = 5.0e-11;
    double solve_tol = 1.0e-12;

    // Solve for the FF-s
    llsq_driver(xml_out, sys, qsq_val, cov_tol, solve_tol);
  }
  catch(std::bad_cast) 
  {
    cerr << argv[0] << ": caught cast error" << endl;
    exit(1);
  }
  catch(const std::string& e) 
  {
    cerr << argv[0] << ": Caught Exception: " << e << endl;
    exit(1);
  }
  catch(std::exception& e) 
  {
    cerr << argv[0] << ": Caught standard library exception: " << e.what() << endl;
    exit(1);
  }
  catch(...)
  {
    cerr << argv[0] << ": caught generic exception during measurement" << endl;
    exit(1);
  }

  pop(xml_out);  // MesonFF
  xml_out.close();
}
