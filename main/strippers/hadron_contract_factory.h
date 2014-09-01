// -*- C++ -*-
// $Id: hadron_contract_factory.h,v 2.0 2008/12/05 04:44:03 edwards Exp $
/*! \file
 *  \brief Factory for producing hadron correlator objects
 */

#ifndef __hadron_2pt_factory_h__
#define __hadron_2pt_factory_h__

#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "io/adat_xmlio.h"
#include "hadron_contract.h"

namespace Strippers
{
  // namespace composition
  using namespace Util;
  using namespace ADATXML;

  //! Hadron 2pt factory (foundry)
  /*! @ingroup hadron */
  typedef SingletonHolder< 
    ObjectFactory<HadronContract, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, int nbin),
		  HadronContract* (*)(XMLReader&,
				      const std::string&,
				      int nbin),
		  StringFactoryError> >
  TheHadronContractFactory;

}


#endif
