// -*- C++ -*-
// $Id: displace_funcmap.h,v 2.0 2008/12/05 04:44:03 edwards Exp $
/*! \file
 *  \brief Displacement function map
 */

#ifndef __displace_funcmap_h__
#define __displace_funcmap_h__

#include "adat/singleton.h"
#include "adat/funcmap.h"
#include "strippers.h"

namespace Strippers
{
  using namespace Util;

  struct DisplacementNames_t
  {
    std::string  particle;
    std::string  wavetype;
  };


  //! Displacement function map 
  typedef SingletonHolder< 
    FunctionMap<DisplacementNames_t, 
		std::string,
		TYPELIST_1(const std::string&),
		DisplacementNames_t (*)(const std::string&),
		StringFunctionMapError> >
  TheDisplacementFuncMap;

  namespace DisplacementCallMapEnv
  { 
    extern bool registered;   // forward decl
  }

} // end namespace Strippers


#endif
