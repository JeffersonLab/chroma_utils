// $Id: build_cache.cc,v 2.0 2008/12/05 04:43:47 edwards Exp $
// Build up a cache file from a user defined list of inputs
// For now, this code is geared towards a Building-Blocks 3pt cache manager

#include "io/adat_xmlio.h"
#include "io/adat_xml_group_reader.h"
#include "adat/handle.h"
#include "ensem/ensem.h"
#include "formfac/formfac_manage_factory.h"
#include "formfac/formfac_manage_aggregate.h"
#include <list>
#include <iostream>

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

	if (arg.quark < 0 || arg.quark >= 2)
	{
	  std::cerr << __func__ << ": snk out of bounds\n";
	  exit(1);
	}
	char ud[2] = {'U', 'D'};

	char filen[1024];
	sprintf(filen, pattern.c_str(), 
		cfg,
		ud[arg.quark-1],
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

  //! 3-pt object
  ADAT::Handle<Manage3PtFuncMap> threept;

  //! User keys
  Array<ThreePtArg> keys;

  try
  {
    XMLReader xml_in(argv[1]);
    XMLReader xml_build(xml_in, record_root);

    read(xml_build, "Keys", keys);

    string threept_id;
    read(xml_build, "ThreePt/ThreePtId", threept_id);
    
    threept = TheManage3PtFuncMapFactory::Instance().createObject(
      threept_id,
      xml_build,
      "ThreePt");
  }
  catch(const std::string& e) 
  {
    cerr << "Caught Exception reading input XML: " << e << endl;
    exit(1);
  }

  cout << argv[0] << ": Loop over keys and build cache file" << endl;

  for(int k=0; k < keys.size(); ++k)
  {
    cout << "Build key " << k << endl;
    EnsemVectorComplexF three = (*threept)[keys[k]];
  }

  // Test the remove
  std::cout << "Erase from map" << std::endl;
//  int cnt = threept.eraseUnused();
//  std::cout << "Erase from map cnt=" << cnt << "\n";

  cout << argv[0] << ": will now exit - look for a cache file" << endl;
}
