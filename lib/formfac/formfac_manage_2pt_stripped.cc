// $Id: formfac_manage_2pt_stripped.cc,v 2.0 2008/12/05 04:43:35 edwards Exp $
/*! \file
 * \brief Manage 2-pt funcs using stripped ensemble files
 */

#include "formfac/formfac_manage_2pt_stripped.h"
#include "formfac/formfac_qsq.h"
#include "formfac/formfac_manage_factory.h"

namespace FF
{

  //----------------------------------------------------------------------------------
  // Parameters
  Manage2PtFuncStrippedParams_t::Manage2PtFuncStrippedParams_t(XMLReader& xml, const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "avg_2pt_func", avg_2pt_func);
    read(paramtop, "LatticeParam", lattice);
    state  = ADATXML::readXMLGroup(paramtop, "State", "StateType");
  }


  // Parameters
  void Manage2PtFuncStrippedParams_t::writeXML(XMLWriter& xml, const std::string& path) const
  {
    push(xml, path);

    write(xml, "avg_2pt_func", avg_2pt_func);
    write(xml, "LatticeParam", lattice);
    xml << state.xml;

    pop(xml);
  }


  // Parameters
  void read(XMLReader& xml, const std::string& path, Manage2PtFuncStrippedParams_t& param)
  {
    Manage2PtFuncStrippedParams_t tmp(xml, path);
    param = tmp;
  }


  // Parameters
  void write(XMLWriter& xml, const std::string& path, const Manage2PtFuncStrippedParams_t& param)
  {
    param.writeXML(xml, path);
  }


  //----------------------------------------------------------------------------------
  // Constructor
  Manage2PtFuncStripped::Manage2PtFuncStripped(ADAT::Handle<State2PtFunc> state_, 
					       bool avg_2pt_func_,
					       const LatticeParam& lattice)
  {
    create(state_, avg_2pt_func_, lattice);
  }


  // Constructor
  void Manage2PtFuncStripped::create(ADAT::Handle<State2PtFunc> state_, 
				     bool avg_2pt_func_,
				     const LatticeParam& lattice)
  {
    state = state_;
    avg_2pt_func = avg_2pt_func_;
    ensemble_info.nbin = 0;
    ensemble_info.lattice = lattice;
  }


  // Read files
  void Manage2PtFuncStripped::do_read2pt(const TwoPtArg& arg)
  {
    //  cout << __func__ << std::endl;

    // Enter if not found
    if (twopt.find(arg) == twopt.end())
    {
      std::string filename = (*state)(0,arg);
      EnsemVectorReal two;
    
      read(filename, two);
      twopt.insert(std::make_pair(arg,two));
    }
  }


  //! Read files
  EnsemVectorReal Manage2PtFuncStripped::operator[](const TwoPtArg& arg)
  {
    TwoPtArg ar(arg);
    if (avg_2pt_func)
      ar.mom = canonicalOrder(arg.mom);

    // Fetch if not found
    if (twopt.find(ar) == twopt.end())
      do_read2pt(ar);

    return twopt.find(ar)->second;
  }


  // Insert a key,value
//  void Manage2PtFuncStripped::insert(const TwoPtArg& arg, const EnsemVectorReal& v)
//  {
//    twopt.insert(std::make_pair(arg,v));
//  }


  // Does key exist?
  bool Manage2PtFuncStripped::exist(const TwoPtArg& arg)
  {
    return (twopt.find(arg) == twopt.end()) ? false : true;
  }


  //! Manager factory
  namespace FormfacManage2PtFuncStrippedEnv
  { 
    //! Anonymous namespace
    namespace
    {
      bool registered = false;

      const std::string name("ENERGY_STRIPPED");

      //-------------------- callback functions ---------------------------------------

      //! Callback
      Manage2PtFuncMap* createFuncMap(XMLReader& xml_in,
				      const std::string& path)
      {
	Manage2PtFuncStrippedParams_t params(xml_in, path);

	// Construct a state
	ADAT::Handle<State2PtFunc> state;
	{
	  std::istringstream  xml_s(params.state.xml);
	  XMLReader  optop(xml_s);
	
	  state = TheState2PtFuncFactory::Instance().createObject(params.state.id,
								  optop,
								  params.state.path);
	}

	return new Manage2PtFuncStripped(state, 
					 params.avg_2pt_func,
					 params.lattice);
      }


      //! Callback
      Manage2PtFunc* createFunc(XMLReader& xml_in,
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
	success &= TheManage2PtFuncFactory::Instance().registerObject(name, createFunc);
	success &= TheManage2PtFuncMapFactory::Instance().registerObject(name, createFuncMap);
	
	registered = true;
      }

      return success;
    }

  }  // end namespace FormfacManageEnv


} // namespace FF
