// $Id: hadron_contract.cc,v 2.0 2008/12/05 04:44:03 edwards Exp $
/*! \file
 *  \brief read hadron contractions
 */

#include "hadron_contract.h"

namespace Strippers
{

#if 0
// NOT SURE IF ANY FUNCTIONS ARE NEEDED HERE YET. THESE WOULD BE
// CATEGORY DEFAULTS

  // Get source location
  // Default version
  Array<int>
  HadronContract::getTSrce(const Array<ForwardProp_t>& forward_headers) const
  {
    Array<int> t_srce = forward_headers[0].source_header.getTSrce();
    int j_decay = forward_headers[0].source_header.j_decay;

    for(int loop=1; loop < forward_headers.size(); ++loop)
    {
      Array<int> t_srce_b = forward_headers[loop].source_header.getTSrce();

      if (t_srce != t_srce_b)
      {
	QDPIO::cerr << __func__ << ": the t_srce in the forward props are not all equal"
		    << endl;
	QDP_abort(1);
      }
    }

    return t_srce;
  }


  // Get source location
  // Default version
  int
  HadronContract::getDecayDir(const Array<ForwardProp_t>& forward_headers) const
  {
    int j_decay = forward_headers[0].source_header.j_decay;

    for(int loop=1; loop < forward_headers.size(); ++loop)
    {
      int j_decay_b = forward_headers[loop].source_header.j_decay;

      if (j_decay != j_decay_b)
      {
	QDPIO::cerr << __func__ << ": the decay_dir in the forward props are not all equal"
		    << endl;
	QDP_abort(1);
      }
    }

    return j_decay;
  }


  //
  // Extract quark prop headers
  // Default version
  //
  ForwardProp_t
  HadronContract::readForwardPropHeader(const std::string& prop_id) const
  {
    ForwardProp_t header;

    try
    {
// No snarfing 
//	ret.prop = TheNamedObjMap::Instance().getData<T>(prop_id);
	
      // Snarf the source info. This is will throw if the source_id is not there
      XMLReader prop_file_xml, prop_record_xml;
      TheNamedObjMap::Instance().get(prop_id).getFileXML(prop_file_xml);
      TheNamedObjMap::Instance().get(prop_id).getRecordXML(prop_record_xml);
   
      // Try to invert this record XML into a ChromaProp struct
      read(prop_record_xml, "/SinkSmear", header);
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << __func__ << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << __func__ << ": caught error = " << e 
		  << endl;
      QDP_abort(1);
    }

    return header;
  }

#endif


} // end namespace Strippers
