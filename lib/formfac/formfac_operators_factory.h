// -*- C++ -*-
// $Id: formfac_operators_factory.h,v 2.0 2008/12/05 04:43:36 edwards Exp $
/*! \file
 *  \brief Factory for producing particle operators
 */

#ifndef __operators_factory_h__
#define __operators_factory_h__

#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "operators.h"

namespace FF
{
  using namespace Util;

  //! Meson operator factory (foundry)
  typedef SingletonHolder< 
    ObjectFactory<Operator<Complex>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  Operator<Complex>* (*)(XMLReader&,
					 const std::string&), StringFactoryError> >
  TheMesonOperatorFactory;


  namespace MesonOperatorEnv
  { 
    extern const bool registered;   // forward decl
  }

}

#endif
