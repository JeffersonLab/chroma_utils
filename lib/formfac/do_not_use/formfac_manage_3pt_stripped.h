// -*- C++ -*-
// $Id: formfac_manage_3pt_stripped.h,v 2.0 2008/12/05 04:43:36 edwards Exp $

/*! \file
 * \brief Manage 3-pt funcs
 */

#ifndef __formfac_manage_3pt_stripped_h__
#define __formfac_manage_3pt_stripped_h__

#include "formfac/formfac_manage_3pt_cache.h"
#include "adat/handle.h"

namespace FF
{

  //----------------------------------------------------------------------------
  //! Class to hold 3-pt functions specific to use individual stripped files
  class Manage3PtFuncStripped : public Manage3PtFuncCache
  {
  public:
    //! Constructor
    Manage3PtFuncStripped(ADAT::Handle<State3PtFunc> state_, 
			  const std::string& cache_file_,
			  int max_map_mb_ = 0);

    //! Destructor
    ~Manage3PtFuncStripped();

  protected:
    //! Read stripped 3pt files
    void do_read3pt(const ThreePtArg& arg);

  private:
    ADAT::Handle<State3PtFunc>  state;
  };

} // namespace FF

#endif
