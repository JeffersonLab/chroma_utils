// -*- C++ -*-
// $Id: hadron_contract.h,v 2.0 2008/12/05 04:44:03 edwards Exp $
/*! \file
 *  \brief Read in hadron correlators
 */

#ifndef __hadron_contract_h__
#define __hadron_contract_h__

#include "io/adat_io.h"
#include "io/adat_xmlio.h"
#include "adat/handle.h"
#include <list>

namespace Strippers
{
  // namespace composition
  using namespace ADAT;
  using namespace ADATIO;
  using namespace ADATXML;
  using XMLArray::Array;

  //! The result of hadron contractions
  /*! @ingroup hadron */
  struct HadronContractResult_t
  {
    std::string         xml;    /*!< XML about each corr group - used to drive the stripper */
    BinaryBufferReader  bin;    /*!< Holds momentum projected correlators */
  };
  

  //! Read hadron correlators
  /*! @ingroup hadron
   *
   * Supports reading of hadron correlators
   */
  class HadronContract
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~HadronContract() {}

    //! Read the correlators
    virtual void operator()(Handle<HadronContractResult_t> had_cont, int ibin) = 0;
    
  protected:
    //! Convenience function to read propagator
//    virtual ForwardProp_t readForwardPropHeader(const std::string& prop_id) const;

    //! Convenience function to get t_srce from headers
//    virtual multi1d<int> getTSrce(const multi1d<ForwardProp_t>& forward_headers) const;

    //! Convenience function to get decay_dir from headers
//    virtual int getDecayDir(const multi1d<ForwardProp_t>& forward_headers) const;
  };

}


#endif
