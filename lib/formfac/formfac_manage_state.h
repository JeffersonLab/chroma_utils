// -*- C++ -*-
// $Id: formfac_manage_state.h,v 2.0 2008/12/05 04:43:36 edwards Exp $
/*! \file
 * \brief Base class for state functions used in manage objects
 */

#ifndef __formfac_manage_state_h__
#define __formfac_manage_state_h__

#include <string>

namespace FF
{
  //----------------------------------------------------------------------------
  //! State function
  template<typename T>
  class StateFunc
  {
  public:
    //! Help with cleanup
    virtual ~StateFunc() {}

    //! Return 
    virtual std::string operator()(int cfg, const T& arg) const = 0;
  };


} // namespace FF

#endif
