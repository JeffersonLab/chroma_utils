// $Id: formfac_manage_3pt.cc,v 2.0 2008/12/05 04:43:35 edwards Exp $
/*! \file
 * \brief Manage 3-pt funcs
 */

#include "formfac/formfac_manage_3pt.h"

namespace FF
{
  //----------------------------------------------------------------------------
  // Return a "view" to a managed object using a prototype key "arg"
  Manage3PtFuncReduced* Manage3PtFunc::createView(const ThreePtArg& prototype_key)
  {
    return new Manage3PtFuncReduced(this, prototype_key);
  }


  //----------------------------------------------------------------------------------
  // Constructor
  Manage3PtFuncReduced::Manage3PtFuncReduced(ADAT::Handle<Manage3PtFunc> obj_, 
					     const ThreePtArg& proto_key_) :
    obj(obj_), proto_key(proto_key_)
  {}

  // Destructor
  Manage3PtFuncReduced::~Manage3PtFuncReduced() {}

  // Return the prototype key
  ThreePtArg Manage3PtFuncReduced::getProtoKey() const {return proto_key;}

  // Return a value given a key
  EnsemVectorComplexF 
  Manage3PtFuncReduced::Manage3PtFuncReduced::operator[](const ThreePtArgReduced& arg)
  {
    // Convert the input argument key with the prototype key
    ThreePtArg key;
    mergeKeys(key, arg, proto_key);

    return (*obj)[key];
  }

  // Number of bins
  int Manage3PtFuncReduced::size() const
  {
    return obj->size();
  }

  // Time extent
  int Manage3PtFuncReduced::timeLen() const
  {
    return obj->timeLen();
  }

  // Decay direction
  int Manage3PtFuncReduced::decayDir() const
  {
    return obj->decayDir();
  }

  // Lattice size
  const ArrayInt& Manage3PtFuncReduced::lattSize() const
  {
    return obj->lattSize();
  }

} // namespace FF
