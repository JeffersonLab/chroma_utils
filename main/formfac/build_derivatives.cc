// $Id: build_derivatives.cc,v 2.0 2008/12/05 04:43:47 edwards Exp $
//
// Build up derivatives
// For now, this code is geared towards a Building-Blocks 3pt cache manager

#include "io/adat_xmlio.h"
#include "io/adat_xml_group_reader.h"
#include "adat/handle.h"
#include "ensem/ensem.h"
#include "formfac/formfac_manage_factory.h"
#include "formfac/formfac_manage_aggregate.h"
#include "formfac/formfac_manage_3pt_cache.h"
#include "parton/parton_derivatives.h"
#include <list>
#include <iostream>

extern "C"
{
#include <lime.h>
}


using namespace ENSEM;
using namespace FF;
using namespace std;


//! Anonymous namespace
namespace
{
  std::ostream& operator<<(std::ostream& s, const Array<int>& d)
  {
    if (d.size() > 0)
    {
      s << d[0];
      for(int i=1; i < d.size(); ++i)
	s << " " << d[i];
    }

    return s;
  }


  std::ostream& operator<<(std::ostream& s, const Array<double>& d)
  {
    if (d.size() > 0)
    {
      s << d[0];
      for(int i=1; i < d.size(); ++i)
	s << " " << d[i];
    }

    return s;
  }
} // namespace



//using namespace ADATXML;
using namespace ENSEM;
using namespace FF;
using namespace std;

// Linkage hack
bool linkage()
{
  bool success = true;
  // Register all factories
  success &= FormfacManage3PtFuncEnv::registerAll();

  return success;
}


int main(int argc, char *argv[])
{
  // Register all factories
  linkage();
  
  // Grab arguments
  if (argc != 2)
  {
    cerr << "Usage:"
	 << argv[0] << ":  <input xml file>" << endl;
    exit(1);
  }

  string record_root = "/BuildCache";

  // Structure of GammaList is: \gamma[1..4], -\gamma5 \gamma[1..4], I \sigma[1,2]...I \sigma[3,4]
  // \sigma only valid for \mu \neq \nu !!!!
  // The signs multiply the gammaList matrix elements
  Array<int> gammaList(14);
  gammaList[ 0]=  1; gammaList[ 1]=  2; gammaList[ 2]=  4; gammaList[ 3]=  8;
  gammaList[ 4]= 14; gammaList[ 5]= 13; gammaList[ 6]= 11; gammaList[ 7]=  7;
  gammaList[ 8]=  3; gammaList[ 9]=  5; gammaList[10]=  9; 
  gammaList[11]=  6; gammaList[12]= 10; gammaList[13]= 12;

  Array<int> gammaSignList(14);
  gammaSignList[ 0]=  1; gammaSignList[ 1]=  1; gammaSignList[ 2]=  1; gammaSignList[ 3]=  1;
  gammaSignList[ 4]= -1; gammaSignList[ 5]=  1; gammaSignList[ 6]= -1; gammaSignList[ 7]=  1;
  gammaSignList[ 8]= -1; gammaSignList[ 9]= -1; gammaSignList[10]= -1; 
  gammaSignList[11]= -1; gammaSignList[12]= -1; gammaSignList[13]= -1;

  //! 3-pt object
  ADAT::Handle<Manage3PtFuncReduced> threept;

  //! User keys
  Array<PiPf> pi_pf;

  //! Ensemble info
  EnsembleInfo ensem_info;

  //! Derivative output
  string output_file;

  try
  {
    XMLReader xml_in(argv[1]);
    XMLReader xml_build(xml_in, record_root);

    read(xml_build, "PiPf", pi_pf);

    string threept_id;
    read(xml_build, "ThreePt/ThreePtId", threept_id);
    
    threept = TheManage3PtFuncReducedFactory::Instance().createObject(
      threept_id,
      xml_build,
      "ThreePt");

    // This should be in there
    read(xml_build, "ThreePt/LatticeParam", ensem_info.lattice);
    ensem_info.nbin = threept->size();

    // Derivative output
    read(xml_build, "output_file", output_file);
  }
  catch(const std::string& e) 
  {
    cerr << "Caught Exception reading input XML: " << e << endl;
    exit(1);
  }

  cout << argv[0] << ": Loop over pi_pf and build cache file" << endl;

  // Use a lime output
  // Begin by opening the file
  FILE *fp;
  if((fp = fopen(output_file.c_str(), "wb")) == NULL)
  {
    std::cerr << __func__ << ": failed to open for writing cache file " 
	      << output_file << std::endl;
    exit(1);
  }

  // Now allocate a lime writer to write out the data
  LimeWriter* writer = limeCreateWriter(fp);
  if( writer == (LimeWriter *)NULL )  
  {
    std::cerr << __func__ << ": unable to open LimeWriter " << std::endl;
    exit(1);
  }

  // Define a bucket to contain the XML data. 
  
  {
    XMLBufferWriter xml;
    write(xml, "Ensem", ensem_info);

    write_cache_entry(writer, 1, 0, xml, "Cache File Header");
  }


  for(int k=0; k < pi_pf.size(); ++k)
  {
    cout << "pi= " << pi_pf[k].p_i << "  pf= " << pi_pf[k].p_f << endl;

    for(int igamma=0; igamma < gammaList.size(); ++igamma)
    {
      for(int mu1=0; mu1 < 4; ++mu1)
      {
	for(int mu2=0; mu2 < 4; ++mu2)
	{
	  cout << "ig= " << igamma << "  mu1= " << mu1 << endl;

	  EnsemVectorComplexF q = 
	    Real(gammaSignList[igamma])*QBar_Gamma_DD_Q(*threept, pi_pf[k], gammaList[igamma], mu1, mu2);

	  // Write a record
	  {
	    ThreePtArgReduced arg;

	    arg.pi_pf = pi_pf[k];
	    arg.gamma = igamma;
	    arg.links.resize(2);
	    arg.links[0] = mu1;
	    arg.links[1] = mu2;

	    ThreePtArg key;
	    mergeKeys(key, arg, threept->getProtoKey());

	    XMLBufferWriter xml;
	    write(xml, "ThreePtArg", key);
	    write_cache_entry(writer, 0, 0, xml, "XML Header");
	    write_cache_entry(writer, 0, 0, q, "Ensemble");
	  }
	}
      }
    }
  }

  limeDestroyWriter(writer);
  fclose(fp);

  cout << argv[0] << ": will now exit - look for a cache file" << endl;
}
