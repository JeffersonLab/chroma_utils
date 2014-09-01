// -*- C++ -*-
// $Id: formfac_manage_3pt_bb.h,v 2.0 2008/12/05 04:43:35 edwards Exp $
/*! \file
 * \brief Manage 3-pt funcs
 */

#ifndef __formfac_manage_3pt_bb_h__
#define __formfac_manage_3pt_bb_h__

#include "formfac/formfac_manage_3pt_cache.h"
#include "formfac/formfac_manage_factory.h"
#include "adat/handle.h"
#include "io/adat_xml_group_reader.h"

namespace FF
{

  //! Manager factory
  namespace FormfacManage3PtFuncCacheBBEnv
  { 
    // Register the callbacks
    bool registerAll(void);
  }


  //----------------------------------------------------------------------------
  //! BB parameters
  struct Manage3PtFuncCacheBBParams_t
  {
    Manage3PtFuncCacheBBParams_t(XMLReader& xml_in, const std::string& path);
    void writeXML(XMLWriter& xml_in, const std::string& path) const;

    GroupXML_t     state;            /*!< xml holding state group */
    std::string    cache_file;
    std::string    cfg_file;
    int            max_map_mb;
    LatticeParam   lattice;
  };

  // Parameters
  void read(XMLReader& xml_in, const std::string& path, Manage3PtFuncCacheBBParams_t& param);

  // Parameters
  void write(XMLWriter& xml_out, const std::string& path, const Manage3PtFuncCacheBBParams_t& param);


  //----------------------------------------------------------------------------
  //! Class to hold 3-pt functions specific to building blocks
  class Manage3PtFuncCacheBB : public Manage3PtFuncCache
  {
  public:
    //! Constructor
    Manage3PtFuncCacheBB(ADAT::Handle<State3PtFunc> state_, 
			 const std::string& cache_file_, 
			 const std::string& cfg_file_,
			 const LatticeParam& lattice,
			 int max_map_mb_ = 0);

    //! Destructor
    ~Manage3PtFuncCacheBB();

  protected:
    //! Read building-blocks 3pt files
    void do_read3pt(const ThreePtArg& arg);

  private:
    ADAT::Handle<State3PtFunc>  state;
    std::list<int>        cfg_list;
  };

} // namespace FF

#endif
