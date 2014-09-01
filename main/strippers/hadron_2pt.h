// -*- C++ -*-
// $Id: hadron_2pt.h,v 2.0 2008/12/05 04:44:03 edwards Exp $
/*! \file
 *  \brief Read hadron 2pt correlators
 */

#ifndef __hadron_2pt_h__
#define __hadron_2pt_h__

#include "ensem/ensem.h"
#include "hadron_contract.h"
//#include "util/ft/sftmom.h"

namespace Strippers
{
  // namespace composition
  using namespace ENSEM;

  //! Read hadron 2pt correlators
  /*! @ingroup hadron */
  struct Hadron2PtCorrs_t
  {
    //! Momentum projected correlator
    struct Mom_t
    {
      Array<int>       mom;    /*!< D-1 momentum of this correlator*/
      Array<Complex>   corr;   /*!< Momentum projected correlator */
    };

    std::string         xml;    /*!< XML about each corr group - used to drive the stripper */
    std::list<Mom_t>    corrs;  /*!< Holds momentum projected correlators */
  };
  

  //! Read hadron 2pt correlators
  /*! @ingroup hadron
   *
   * Supports reading of hadron 2pt correlators
   *
   */
  class Hadron2PtCorr : public HadronContract
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~Hadron2PtCorr() {}

    //! Read the correlators
    virtual void operator()(Handle<HadronContractResult_t> had_cont, int ibin) = 0;
    
  protected:
    //! Take a record and convert it back to a correlator object
    virtual Handle<Hadron2PtCorrs_t> unserialize(Handle<HadronContractResult_t> had_cont) const;
  };

}


#endif
