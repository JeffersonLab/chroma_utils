// $Id: hadron_contract_aggregate.cc,v 2.0 2008/12/05 04:44:03 edwards Exp $
/*! \file
 *  \brief All hadron contraction constructors
 */

#include "hadron_contract_aggregate.h"

#include "simple_meson_2pt_w.h"
//#include "simple_baryon_2pt_w.h"

namespace Strippers
{

  //! Registration aggregator
  namespace HadronContractEnv
  {
    namespace
    {
      //! Local registration flag
      bool registered = false;
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Hadron
	success &= SimpleMeson2PtEnv::registerAll();
//	success &= SimpleBaryon2PtEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
