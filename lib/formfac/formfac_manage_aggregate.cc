// $Id: formfac_manage_aggregate.cc,v 2.0 2008/12/05 04:43:36 edwards Exp $
/*! \file
 *  \brief All formfactor factories
 */

#include "formfac/formfac_manage_aggregate.h"

#include "formfac/formfac_manage_2pt_stripped.h"
#include "formfac/formfac_manage_E_stripped.h"
#include "formfac/formfac_manage_Z_stripped.h"

#include "formfac/formfac_manage_3pt_db.h"

namespace FF
{

  //! Registration aggregator
  namespace FormfacManage2PtFuncEnv
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
	success &= FormfacManage2PtFuncStrippedEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }


  //! Registration aggregator
  namespace FormfacManageEnergyFuncEnv
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
	success &= FormfacManageEnergyFuncStrippedEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }


  //! Registration aggregator
  namespace FormfacManageAmpFuncEnv
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
	success &= FormfacManageAmpFuncStrippedEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }


  //! Registration aggregator
  namespace FormfacManage3PtFuncEnv
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
	success &= FormfacManage3PtFuncDBEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
