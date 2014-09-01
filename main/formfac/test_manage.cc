// $Id: test_manage.cc,v 2.0 2008/12/05 04:43:49 edwards Exp $
// Vector-pseudoscalar form-factor code using fitting method

#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "formfac/formfac.h"
#include "state.h"
#include <list>
#include <iostream>

#include <sys/stat.h>
#include <sys/types.h>

namespace FF
{
  //! Make a state
  class RhoPion3PtFunc : public State3PtFunc
  {
  public:
    //! Constructor
    RhoPion3PtFunc(const std::string& pattern_) : pattern(pattern_) {}

    //! Construct 3pt names
    std::string operator()(int cfg, const ThreePtArg& arg) const
      {
	ArrayInt p_i = arg.pi_pf.p_i;
	ArrayInt p_f = arg.pi_pf.p_f;
	ArrayInt q   = p_f - p_i;
  
	// cout << __func__ << ": pattern=" << pattern.c_str() << std::endl;

	char qxsign = (q[0] < 0 ) ? '-' : '+';
	char qysign = (q[1] < 0 ) ? '-' : '+';
	char qzsign = (q[2] < 0 ) ? '-' : '+';

	char pfxsign = (p_f[0] < 0 ) ? '-' : '+';
	char pfysign = (p_f[1] < 0 ) ? '-' : '+';
	char pfzsign = (p_f[2] < 0 ) ? '-' : '+';

	char xyz[3] = {'X', 'Y', 'Z'};

	char filen[1024];
	sprintf(filen, pattern.c_str(),
		pfzsign, abs(p_f[2]), pfysign, abs(p_f[1]), pfxsign, abs(p_f[0]),
		cfg,
		xyz[arg.src[0]],
		qzsign, abs(q[2]), qysign, abs(q[1]), qxsign, abs(q[0]));

	std::cout << __func__ << ": filen=" << filen << "\n";

	return std::string(filen);
      }

  private:
    std::string pattern;
  };

} // namespace FF



//using namespace ADATXML;
using namespace ENSEM;
using namespace FF;
using namespace std;


int main(int argc, char *argv[])
{
  string pattern = "/cache/HASTE/NF0/iso/6p5_24_48_wl/bb.unpack/cl/data.5p0_5p0.dir/pz%c%1d_py%c%1d_px%c%1d/%d/RHO_%c_1-PION_1_%c%1d_%c%1d_%c%1d.bb";
  string cache_file = "my_cache.lime";
  string cfg_file = "/cache/HASTE/NF0/iso/6p5_24_48_wl/bb.unpack/cl/data.5p0_5p0.dir/cfg_list";

  ADAT::Handle<State3PtFunc> state(new RhoPion3PtFunc(pattern));
  int max_map_mb = 1;
  Manage3PtFuncBB threept(state, cache_file, cfg_file, max_map_mb);

  ArrayInt one(3);
  one[0] = 1;
  one[1] = 0;
  one[2] = 0;

  for(int p=0; p < 3; ++p)
  {
    for(int g=0; g < 15; ++g)
    {
      ThreePtArg  arg;
      arg.pi_pf.p_i = one;
      arg.pi_pf.p_f = one;
      arg.g      = g;
      arg.src.resize(1);
      arg.src[0] = p;

      EnsemVectorComplexF& three = threept[arg];

      if (g == 8 && p == 0)
	write("my_3pt.dat", three);
    }
  }

  // Test the remove
  std::cout << "Erase from map" << std::endl;
  int cnt = threept.eraseUnused();
  std::cout << "Erase from map cnt=" << cnt << "\n";
}
