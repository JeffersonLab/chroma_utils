// -*- C++ -*-
// $Id: current_ops_factory.h,v 2.0 2008/12/05 04:43:47 edwards Exp $
/*! \file
 *  \brief Factory for current operators
 */

#ifndef __current_operators_factory_h__
#define __current_operators_factory_h__

#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "current_ops.h"

namespace FF
{
  using namespace Util;

  //! Meson operator factory (foundry)
  typedef SingletonHolder< 
    ObjectFactory<CurrentOperator, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  CurrentOperator* (*)(XMLReader&,
				       const std::string&), StringFactoryError> >
  TheCurrentOperatorFactory;


  namespace CurrentOperatorEnv
  { 
    bool registerAll();
  }

}

#endif
