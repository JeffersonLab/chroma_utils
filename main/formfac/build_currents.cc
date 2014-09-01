// $Id: build_currents.cc,v 2.0 2008/12/05 04:43:47 edwards Exp $
//
// Build up derivatives
// For now, this code is geared towards a Building-Blocks 3pt cache manager

#include "io/adat_xmlio.h"
#include "io/adat_xml_group_reader.h"
#include "adat/handle.h"
#include "ensem/ensem.h"
#include "formfac/formfac.h"
#include "formfac/formfac_manage_aggregate.h"
#include "parton/parton_derivatives.h"
#include <list>
#include <iostream>

extern "C"
{
#include <lime.h>
}

namespace FF
{
  //! Make a state appropriate to the DWF/DWF BB's
  class Local3PtFunc : public State3PtFunc
  {
  public:
    //! Constructor
    Local3PtFunc(const std::string& pattern_) : pattern(pattern_) {}

    //! Construct 3pt names
    //! Construct 3pt names
    std::string operator()(int cfg, const ThreePtArg& arg) const
      {
	ArrayInt p_i = arg.pi_pf.p_i;
	ArrayInt p_f = arg.pi_pf.p_f;
	ArrayInt q   = p_f - p_i;      // note sign convention

//	std::cout << __func__ << ": pattern=" << pattern.c_str() << std::endl;

	char qxsign = (q[0] < 0 ) ? '-' : '+';
	char qysign = (q[1] < 0 ) ? '-' : '+';
	char qzsign = (q[2] < 0 ) ? '-' : '+';

	char pfxsign = (p_f[0] < 0 ) ? '-' : '+';
	char pfysign = (p_f[1] < 0 ) ? '-' : '+';
	char pfzsign = (p_f[2] < 0 ) ? '-' : '+';

	if (arg.snk[0] < 0 || arg.snk[0] >= 2)
	{
	  std::cerr << __func__ << ": snk out of bounds\n";
	  exit(1);
	}

	char filen[1024];
	sprintf(filen, pattern.c_str(), 
		cfg,
		pfzsign, abs(p_f[2]), pfysign, abs(p_f[1]), pfxsign, abs(p_f[0]),
		qzsign, abs(q[2]), qysign, abs(q[1]), qxsign, abs(q[0]));

	std::cout << "3pt filen=" << filen << std::endl;

	return std::string(filen);
      }

  private:
    std::string pattern;
  };

 
  //! Manager factory
  namespace LocalState3PtEnv
  { 
    //! Anonymous namespace
    namespace
    {
      bool registered = false;

      //-------------------- callback functions ---------------------------------------

      //! Callback
      State3PtFunc* createFunc(XMLReader& xml_in,
			       const std::string& path)
      {
	XMLReader statexml(xml_in, path);
	std::string pattern;
	read(statexml, "pattern", pattern);

	return new Local3PtFunc(pattern);
      }

    }  // anonymous namespace


    // Register the callbacks
    bool registerAll(void) 
    {
      bool success = true;

      if (! registered)
      {
	//! Register all the factories
	success &= TheState3PtFuncFactory::Instance().registerObject(std::string("LOCAL_BB"),
								     createFunc);

	registered = true;
      }

      return success;
    }

  }  // end namespace LocalState3PtEnv


  //! Anonymous namespace
  namespace
  {
    std::ostream& operator<<(std::ostream& s, const ArrayInt& d)
    {
      if (d.size() > 0)
      {
	s << d[0];
	for(int i=1; i < d.size(); ++i)
	  s << " " << d[i];
      }

      return s;
    }


    std::ostream& operator<<(std::ostream& s, const ArrayDouble& d)
    {
      if (d.size() > 0)
      {
	s << d[0];
	for(int i=1; i < d.size(); ++i)
	  s << " " << d[i];
      }

      return s;
    }
  }

} // namespace FF



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
  success &= LocalState3PtEnv::registerAll();

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
  ArrayInt gammaList(14);
  gammaList[ 0]=  1; gammaList[ 1]=  2; gammaList[ 2]=  4; gammaList[ 3]=  8;
  gammaList[ 4]= 14; gammaList[ 5]= 13; gammaList[ 6]= 11; gammaList[ 7]=  7;
  gammaList[ 8]=  3; gammaList[ 9]=  5; gammaList[10]=  9; 
  gammaList[11]=  6; gammaList[12]= 10; gammaList[13]= 12;

  ArrayInt gammaSignList(14);
  gammaSignList[ 0]=  1; gammaSignList[ 1]=  1; gammaSignList[ 2]=  1; gammaSignList[ 3]=  1;
  gammaSignList[ 4]= -1; gammaSignList[ 5]=  1; gammaSignList[ 6]= -1; gammaSignList[ 7]=  1;
  gammaSignList[ 8]= -1; gammaSignList[ 9]= -1; gammaSignList[10]= -1; 
  gammaSignList[11]= -1; gammaSignList[12]= -1; gammaSignList[13]= -1;

  //! 3-pt object
  ADAT::Handle<Manage3PtFunc> threept;

  //! User keys
  Array<PiPf> pi_pf;

  //! Ensemble info
  EnsembleInfo ensem_info;

  //! Derivative output
  string output_file;
  string prefix;

  try
  {
    XMLReader xml_in(argv[1]);
    XMLReader xml_build(xml_in, record_root);

    read(xml_build, "Threept/prefix", prefix);
    read(xml_build, "PiPf", pi_pf);

    string threept_id;
    read(xml_build, "ThreePt/ThreePtId", threept_id);
    
    threept = TheManage3PtFuncFactory::Instance().createObject(
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
	cout << "ig= " << igamma << "  mu1= " << mu1 << endl;

	EnsemVectorComplexF u = 
	  Real(gammaSignList[igamma])*UBar_Gamma_D_U(*threept, pi_pf[k], gammaList[igamma], 1 << mu1);

	// Write a record
	{
	  std::ostringstream os;
	  os << "cons_" << prefix << "_" << mu << "_pi" << pi_pf[k].p_i << "_pf" << pi_pf[k].p_f;
	  write(os.str(), u);
	}
      }
    }
  }

  limeDestroyWriter(writer);
  fclose(fp);

  cout << argv[0] << ": will now exit - look for a cache file" << endl;
}
