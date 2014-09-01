// -*- C++ -*-
/*! \file
 * \brief Support for unordered maps
 */

#include "map_traits.h"
#include "ffdb_db.h"
#include "ffdb_hash_func.h"

namespace ADAT
{
  //! Hash function call
  size_t hashfunc(const std::string& that)
  {
    return std::hash<std::string>()(that);
//    return (*__ffdb_default_hash)(that.data(), that.size());
  }

  //! Equivalence
  bool hashcmp(const std::string& that1, const std::string& that2)
  {
    return (that1 == that2) ? true : false;
  }

} // namespace ADAT

