// $Id: formfac_manage_Z_stripped.cc,v 2.0 2008/12/05 04:52:53 edwards Exp $
/*! \file
 * \brief Manage amplitudes
 */

#include "formfac/formfac_manage_Z_stripped.h"
#include "formfac/formfac_qsq.h"
#include "formfac/formfac_manage_factory.h"

namespace FF
{

  //----------------------------------------------------------------------------------
  // Parameters
  ManageAmpFuncStrippedParams_t::ManageAmpFuncStrippedParams_t(XMLReader& xml, const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "avg_2pt_func", avg_2pt_func);
    read(paramtop, "LatticeParam", lattice);
    state  = ADATXML::readXMLGroup(paramtop, "State", "StateType");
  }


  // Parameters
  void ManageAmpFuncStrippedParams_t::writeXML(XMLWriter& xml, const std::string& path) const
  {
    push(xml, path);

    write(xml, "avg_2pt_func", avg_2pt_func);
    write(xml, "LatticeParam", lattice);
    xml << state.xml;

    pop(xml);
  }


  // Parameters
  void read(XMLReader& xml, const std::string& path, ManageAmpFuncStrippedParams_t& param)
  {
    ManageAmpFuncStrippedParams_t tmp(xml, path);
    param = tmp;
  }


  // Parameters
  void write(XMLWriter& xml, const std::string& path, const ManageAmpFuncStrippedParams_t& param)
  {
    param.writeXML(xml, path);
  }


  //----------------------------------------------------------------------------------
  // Constructor
  ManageAmpFuncStripped::ManageAmpFuncStripped(ADAT::Handle<StateAmpFunc> state_,
					       bool avg_2pt_func_,
					       const LatticeParam& lattice)
  {
    create(state_, avg_2pt_func_, lattice);
  }


  // Constructor
  void ManageAmpFuncStripped::create(ADAT::Handle<StateAmpFunc> state_,
				     bool avg_2pt_func_,
				     const LatticeParam& lattice)
  {
    state = state_;
    avg_2pt_func = avg_2pt_func_;
    ensemble_info.nbin = 0;
    ensemble_info.lattice = lattice;
  }


  // Read files
  void ManageAmpFuncStripped::do_readZ(const AmpArg& arg)
  {
    //  cout << __func__ << std::endl;

    // Enter if not found
    if (Z.find(arg) == Z.end())
    {
      std::string filename = (*state)(-1,arg);
      EnsemReal two;
    
      read(filename, two);
      Z.insert(std::make_pair(arg,two));
    }
  }


  //! Read files
  EnsemReal ManageAmpFuncStripped::operator[](const AmpArg& arg)
  {
    AmpArg ar(arg);
    if (avg_2pt_func)
      ar.mom = canonicalOrder(arg.mom);

    // Fetch if not found
    if (Z.find(ar) == Z.end())
      do_readZ(ar);

    return Z.find(ar)->second;
  }


  // Insert a key,value
  void ManageAmpFuncStripped::insert(const AmpArg& arg, const EnsemReal& v)
  {
    Z.insert(std::make_pair(arg,v));
  }


  // Does key exist?
  bool ManageAmpFuncStripped::exist(const AmpArg& arg)
  {
    return (Z.find(arg) == Z.end()) ? false : true;
  }


  //! Manager factory
  namespace FormfacManageAmpFuncStrippedEnv
  { 
    //! Anonymous namespace
    namespace
    {
      bool registered = false;

      const std::string name("ENERGY_STRIPPED");

      //-------------------- callback functions ---------------------------------------

      //! Callback
      ManageAmpFuncMap* createFuncMap(XMLReader& xml_in,
				      const std::string& path)
      {
	ManageAmpFuncStrippedParams_t params(xml_in, path);

	// Construct a state
	ADAT::Handle<StateAmpFunc> state;
	{
	  std::istringstream  xml_s(params.state.xml);
	  XMLReader  optop(xml_s);
	
	  state = TheStateAmpFuncFactory::Instance().createObject(params.state.id,
								  optop,
								  params.state.path);
	}

	return new ManageAmpFuncStripped(state, 
					 params.avg_2pt_func,
					 params.lattice);
      }


      //! Callback
      ManageAmpFunc* createFunc(XMLReader& xml_in,
				const std::string& path)
      {
	return createFuncMap(xml_in, path);
      }

    }  // anonymous namespace


    // Register the callbacks
    bool registerAll(void) 
    {
      bool success = true;

      if (! registered)
      {
	//! Register all the factories
	success &= TheManageAmpFuncFactory::Instance().registerObject(name, createFunc);
	success &= TheManageAmpFuncMapFactory::Instance().registerObject(name, createFuncMap);
	
	registered = true;
      }

      return success;
    }

  }  // end namespace FormfacManageEnv

} // namespace FF
