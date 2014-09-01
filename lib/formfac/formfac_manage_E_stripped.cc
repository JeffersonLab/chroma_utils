// $Id: formfac_manage_E_stripped.cc,v 2.0 2008/12/05 04:52:53 edwards Exp $
/*! \file
 * \brief Manage energy funcs
 */

#include "formfac/formfac_manage_E_stripped.h"
#include "formfac/formfac_qsq.h"
#include "formfac/formfac_manage_factory.h"

namespace FF
{

  //----------------------------------------------------------------------------------
  // Parameters
  ManageEnergyFuncStrippedParams_t::ManageEnergyFuncStrippedParams_t(XMLReader& xml, const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "LatticeParam", lattice);
    state  = ADATXML::readXMLGroup(paramtop, "State", "StateType");
  }


  // Parameters
  void ManageEnergyFuncStrippedParams_t::writeXML(XMLWriter& xml, const std::string& path) const
  {
    push(xml, path);

    write(xml, "LatticeParam", lattice);
    xml << state.xml;

    pop(xml);
  }


  // Parameters
  void read(XMLReader& xml, const std::string& path, ManageEnergyFuncStrippedParams_t& param)
  {
    ManageEnergyFuncStrippedParams_t tmp(xml, path);
    param = tmp;
  }


  // Parameters
  void write(XMLWriter& xml, const std::string& path, const ManageEnergyFuncStrippedParams_t& param)
  {
    param.writeXML(xml, path);
  }


  //----------------------------------------------------------------------------------
  // Constructor
  ManageEnergyFuncStripped::ManageEnergyFuncStripped(ADAT::Handle<StateEnergyFunc> state_,
						     const LatticeParam& lattice)
  {
    create(state_, lattice);
  }


  // Constructor
  void ManageEnergyFuncStripped::create(ADAT::Handle<StateEnergyFunc> state_,
					const LatticeParam& lattice)
  {
    state = state_;
    ensemble_info.nbin = 0;
    ensemble_info.lattice = lattice;
  }


  // Read files
  void ManageEnergyFuncStripped::do_readE(const EnergyArg& p)
  {
    //  cout << __func__ << std::endl;

    // Enter if not found
    if (E.find(p) == E.end())
    {
      std::string filename = (*state)(-1,p);
      EnsemReal e;
    
      read(filename, e);
      E.insert(std::make_pair(p,e));
    }
  }


  //! Read files
  EnsemReal ManageEnergyFuncStripped::operator[](const EnergyArg& arg0)
  {
    EnergyArg arg = arg0;
    arg.mom = canonicalOrder(arg0.mom);

    // Fetch if not found
    if (E.find(arg) == E.end())
      do_readE(arg);

    return E.find(arg)->second;
  }


  // Insert a key,value
  void ManageEnergyFuncStripped::insert(const EnergyArg& arg0, const EnsemReal& v)
  {
    EnergyArg arg = arg0;
    arg.mom = canonicalOrder(arg0.mom);
    E.insert(std::make_pair(arg,v));
  }


  // Does key exist?
  bool ManageEnergyFuncStripped::exist(const EnergyArg& arg0)
  {
    EnergyArg arg = arg0;
    arg.mom = canonicalOrder(arg0.mom);

    return (E.find(arg) == E.end()) ? false : true;
  }


  //! Manager factory
  namespace FormfacManageEnergyFuncStrippedEnv
  { 
    //! Anonymous namespace
    namespace
    {
      bool registered = false;

      const std::string name("ENERGY_STRIPPED");

      //-------------------- callback functions ---------------------------------------

      //! Callback
      ManageEnergyFuncMap* createFuncMap(XMLReader& xml_in,
					 const std::string& path)
      {
	ManageEnergyFuncStrippedParams_t params(xml_in, path);

	// Construct a state
	ADAT::Handle<StateEnergyFunc> state;
	{
	  std::istringstream  xml_s(params.state.xml);
	  XMLReader  optop(xml_s);
	
	  state = TheStateEnergyFuncFactory::Instance().createObject(params.state.id,
								     optop,
								     params.state.path);
	}

	return new ManageEnergyFuncStripped(state, 
					    params.lattice);
      }


      //! Callback
      ManageEnergyFunc* createFunc(XMLReader& xml_in,
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
	success &= TheManageEnergyFuncFactory::Instance().registerObject(name, createFunc);
	success &= TheManageEnergyFuncMapFactory::Instance().registerObject(name, createFuncMap);
	
	registered = true;
      }

      return success;
    }

  }  // end namespace FormfacManageEnv

} // namespace FF
