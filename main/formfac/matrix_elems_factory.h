// -*- C++ -*-
// $Id: matrix_elems_factory.h,v 2.0 2008/12/05 04:43:48 edwards Exp $
/*! \file
 *  \brief Factory for producing matrix elements
 */

#ifndef __matrix_elems_factory_h__
#define __matrix_elems_factory_h__

#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "meson_matrix_elems.h"
#include "baryon_matrix_elems.h"

namespace FF
{
  using namespace Util;

  //! MatrixElement factory (foundry)
  typedef SingletonHolder< 
    ObjectFactory<MesonMatrixElement,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  MesonMatrixElement* (*)(XMLReader&,
					  const std::string&), StringFactoryError> >
  TheMesonMatrixElementFactory;


  //! MatrixElement factory (foundry)
  typedef SingletonHolder< 
    ObjectFactory<BaryonMatrixElement,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  BaryonMatrixElement* (*)(XMLReader&,
					   const std::string&), StringFactoryError> >
  TheBaryonMatrixElementFactory;


  namespace MesonMatrixElementEnv
  { 
    bool registerAll();
  }

}

#endif
