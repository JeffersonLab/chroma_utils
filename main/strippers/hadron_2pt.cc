// $Id: hadron_2pt.cc,v 2.0 2008/12/05 04:44:03 edwards Exp $
/*! \file
 *  \brief Read hadron 2pt correlators
 */

#include "hadron_2pt.h"
//#include "util/ft/sftmom.h"

namespace Strippers
{

  // Un-serialize the structure
  Handle<Hadron2PtCorrs_t>
  Hadron2PtCorr::unserialize(Handle<HadronContractResult_t> had_cont) const
  {
    Handle<Hadron2PtCorrs_t> had_2pt(new Hadron2PtCorrs_t());

    had_2pt->xml = had_cont->xml;

    // Read the momentum projected
    while(! had_cont->bin.eof())
    {
      Hadron2PtCorrs_t::Mom_t moms;

      read(had_cont->bin, moms.mom);
      read(had_cont->bin, moms.corr);

      had_2pt->corrs.push_back(moms);
    }

    return had_2pt;
  }

} // end namespace Strippers
