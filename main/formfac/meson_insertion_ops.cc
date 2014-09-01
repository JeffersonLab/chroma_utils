// $Id: meson_insertion_ops.cc,v 2.0 2008/12/05 04:43:48 edwards Exp $
/*! \file
 * \brief Meson insertion operators
 */

#include "meson_insertion_ops.h"
#include "insertion_ops_factory.h"
#include <iostream>

using namespace std;

namespace FF
{
  //---------------------------------------------------------------------------
  // Full constructor
  Pi0VecCurPi0InsertOp::Pi0VecCurPi0InsertOp(const InsertionOperatorParams& p) : params(p)
  {
#if 0
    // Construct 3pt manage object
    try
    {
      std::istringstream  xml_s(params.threept.xml);
      XMLReader  optop(xml_s);
	
      threept = TheManage3PtFuncFactory::Instance().createObject(params.threept.id,
								 optop,
								 params.threept.path);

      // This should be in there
      read(optop, "ThreePt/LatticeParam", ensem_info.lattice);
      ensem_info.nbin = threept->size();
    }
    catch(const std::string& e) 
    {
      cerr << "Caught Exception reading input XML: " << e << endl;
      exit(1);
    }
#endif
  }


 
  //! Meson sources
  namespace MesonInsertionOperatorEnv
  { 
    //! Anonymous namespace
    namespace
    {
      bool registered = false;

      //-------------------- callback functions ---------------------------------------

      //! pi0|V|pi0
      MesonInsertionOperator* matPi0VecCurPi0(XMLReader& xml_in,
					  const std::string& path)
      {
	return new Pi0VecCurPi0InsertOp(InsertionOperatorParams(xml_in, path));
      }

    }  // anonymous namespace


    // Register all the possible matrix elements
    bool registerAll() 
    {
      bool success = true;

      if (! registered)
      {
	//! Register all the new objects
	success &= TheMesonInsertionOperatorFactory::Instance().registerObject(std::string("pi0|VecCur|pi0"),
									       matPi0VecCurPi0);

	registered = true;
      }

      return success;
    }

  }  // end namespace InsertionOperatorEnv

} // namespace FF
