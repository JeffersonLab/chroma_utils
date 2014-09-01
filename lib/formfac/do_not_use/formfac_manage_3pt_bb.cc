// $Id: formfac_manage_3pt_bb.cc,v 2.0 2008/12/05 04:43:35 edwards Exp $
/*! \file
 * \brief Manage 3-pt funcs
 */

#include "formfac/formfac_manage_3pt_bb.h"
#include "formfac/formfac_bb.h"
#include "formfac/formfac_manage_factory.h"

#undef FF_DEBUG

namespace FF
{

  //----------------------------------------------------------------------------------
  // Parameters
  Manage3PtFuncCacheBBParams_t::Manage3PtFuncCacheBBParams_t(XMLReader& xml, const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "cache_file", cache_file);
    read(paramtop, "cfg_file", cfg_file);
    read(paramtop, "max_map_mb", max_map_mb);
    read(paramtop, "LatticeParam", lattice);
    state  = ADATXML::readXMLGroup(paramtop, "State", "StateType");
  }


  // Parameters
  void Manage3PtFuncCacheBBParams_t::writeXML(XMLWriter& xml, const std::string& path) const
  {
    push(xml, path);

    write(xml, "cache_file", cache_file);
    write(xml, "cfg_file", cfg_file);
    write(xml, "max_map_mb", max_map_mb);
    write(xml, "LatticeParam", lattice);
    xml << state.xml;

    pop(xml);
  }


  // Parameters
  void read(XMLReader& xml, const std::string& path, Manage3PtFuncCacheBBParams_t& param)
  {
    Manage3PtFuncCacheBBParams_t tmp(xml, path);
    param = tmp;
  }


  // Parameters
  void write(XMLWriter& xml, const std::string& path, const Manage3PtFuncCacheBBParams_t& param)
  {
    param.writeXML(xml, path);
  }


  //----------------------------------------------------------------------------------
  //! Constructor
  Manage3PtFuncCacheBB::Manage3PtFuncCacheBB(ADAT::Handle<State3PtFunc> state_, 
					     const std::string& cache_file_,
					     const std::string& cfg_file,
					     const LatticeParam& lattice,
					     int max_map_mb_) :
    Manage3PtFuncCache(cache_file_,max_map_mb_), state(state_)
  {
    EnsembleInfo  local_info;  /*!< Held here until the whole ensemble is initialized*/
    local_info.lattice = lattice;

    // Read cfg list
    std::ifstream fp(cfg_file.c_str(), std::ifstream::in);
    if (fp.fail())
    {
      std::cerr << __func__ << ": error opening file=" << cfg_file << std::endl;
      exit(1);
    }

    int ll;
    fp >> ll;
    while (! fp.eof() && ! fp.fail())
    {
      cfg_list.push_back(ll);
      fp >> ll;
    }

    // If the cache manager is initialized, just keep a local copy of the
    // ensemble info. If it is not initialized, then hold a local copy
    // to finally initialize the whole thing.
    if (this->isInitP())
    {
      local_info = this->ensembleInfo();

      if (local_info.nbin != cfg_list.size())
      {
	std::cerr << __func__ << ": inconsistent cfg_list size and initialized size" << std::endl;
	exit(1);
      }
    }
    else
    {
      local_info.nbin = cfg_list.size();
      initEnsemble(local_info);
    }

    std::cout << __func__ << ": nbins=" << this->size() << std::endl;

    fp.close();
  }


  //! Destructor
  Manage3PtFuncCacheBB::~Manage3PtFuncCacheBB() {}


  //! Read files
  void Manage3PtFuncCacheBB::do_read3pt(const ThreePtArg& arg)
  {
    std::cout << __func__ << ": entering" << std::endl;

    StopWatch swatch,snoop;
    swatch.reset();
    swatch.start();
    snoop.reset();
    snoop.start();

#ifdef FF_DEBUG
    std::cout << __func__ << "find:  arg=" << arg;
#endif

    if (! this->isInitP())
    {
      std::cerr << __func__ << ": ensemble not initialized" << std::endl;
      exit(1);
    }

    //! TR1 unordered_map version (a hash table)
    typedef std::tr1::unordered_map<ThreePtArg, 
      EnsemVectorComplexF, 
      UnorderedMapTraits<ThreePtArg>,
      UnorderedMapTraits<ThreePtArg> > MapType_t;

    MapType_t    threept;

    bool first_configP = true;

    int nb=0;
    for(std::list<int>::const_iterator mm = cfg_list.begin();
	mm != cfg_list.end();
	++mm)
    {
      if (nb >= this->size())
      {
	std::cerr << __func__ << ": internal error - number of cfgs inconsistent" << std::endl;
	exit(1);
      }

      int cfg = *mm;
      std::string filename = (*state)(cfg, arg);

      if (first_configP)
        std::cout << __func__ << ": read " << filename << std::endl;

      BuildingBlocks_t bar;
      read(filename, bar);
    
      for(int l=0; l < bar.links.size(); ++l)
      {
	const Array<int>& link_value = bar.links[l].link_value;

	for(int g=0; g < bar.links[l].gamma.size(); ++g)
	{
	  for(int m=0; m < bar.links[l].gamma[g].momenta.size(); ++m)
	  {
	    const Array<int>& inser_mom    = bar.links[l].gamma[g].momenta[m].inser_mom;
	    const Array<ComplexF>& current = bar.links[l].gamma[g].momenta[m].current;

	    // Time extent needed later. Also hold onto ensemble info.
	    int time_extent = current.size();
	    
	    // The single config correlator
	    VectorComplexF pt;
	    pt.resizeObs(time_extent);

	    for(int t=0; t < current.size(); ++t)
	    {
	      pokeObs(pt, current[t], t);
	    }

	    // Key to the 3pt
	    ThreePtArg ar(arg);
	    ar.pi_pf.p_f = arg.pi_pf.p_f;

	    // ar.pi_pf.p_i = arg.pi_pf.p_f - inser_mom;
	    ar.pi_pf.p_i = arg.pi_pf.p_f + inser_mom; //*** changed to deal with the phase convention change ***

	    ar.gamma       = g;
	    ar.links       = link_value;

	    // Create a zero entry the first time around
	    if (first_configP)
	    {
	      EnsemVectorComplexF thr;
	      thr.resize(this->size());
	      thr.resizeObs(time_extent);
	      thr = zero;

	      threept.insert(std::make_pair(ar,thr));
	    }

	    // Insert into the correct cfg
	    pokeEnsem(threept.find(ar)->second, pt, nb);
	  }
	}
      }

      first_configP = false;
      nb++;
    }

    swatch.stop();
    std::cout << __func__ << ": time to read BB = "
	      << swatch.getTimeInSeconds() 
	      << " secs" << std::endl;
    swatch.reset();
    swatch.start();

    // Insert the desired guy into the map first and mark it so it cannot be erased
    {
      MapType_t::iterator mm = threept.find(arg);
      if (mm != threept.end())
      {
	this->insert(mm->first, mm->second);
	this->mark(arg);
      }
      else
      {
	std::cerr << __func__ << ": internal error: entry not found" << std::endl;
	exit(1);
      }
    }
	    
    // Insert into the map
    for(MapType_t::const_iterator mm = threept.begin(); mm != threept.end(); ++mm)
    {
      this->insert(mm->first, mm->second);
    }

    // Grab and throw away the desired 3pt just to mark it used (bump count).
    EnsemVectorComplexF foo = (*this)[arg];

    swatch.stop();
    snoop.stop();

    std::cout << __func__ << ": time to insert = "
	      << swatch.getTimeInSeconds() 
	      << " secs" << std::endl;

    std::cout << __func__ << ": total time = "
	      << snoop.getTimeInSeconds() 
	      << " secs" << std::endl;

  }



  //! Manager factory
  namespace FormfacManage3PtFuncCacheBBEnv
  { 
    //! Anonymous namespace
    namespace
    {
      bool registered = false;

      const std::string name("THREEPT_CACHE_BB");

      //-------------------- callback functions ---------------------------------------

      //! Callback
      Manage3PtFuncMap* createFuncMap(XMLReader& xml_in,
				      const std::string& path)
      {
	Manage3PtFuncCacheBBParams_t params(xml_in, path);

	// Construct a state
	ADAT::Handle<State3PtFunc> state;
	{
	  std::istringstream  xml_s(params.state.xml);
	  XMLReader  optop(xml_s);
	
	  state = TheState3PtFuncFactory::Instance().createObject(params.state.id,
								  optop,
								  params.state.path);
	}

	return new Manage3PtFuncCacheBB(state, 
					params.cache_file, 
					params.cfg_file, 
					params.lattice, 
					params.max_map_mb);
      }


      //! Callback
      Manage3PtFunc* createFunc(XMLReader& xml_in,
				const std::string& path)
      {
	return createFuncMap(xml_in, path);
      }

    }  // anonymous namespace


    // Register the callbacks
    bool registerAll(void) 
    {
      bool success = true;

      if (! registered)
      {
	//! Register all the factories
	success &= TheManage3PtFuncFactory::Instance().registerObject(name, createFunc);
	success &= TheManage3PtFuncMapFactory::Instance().registerObject(name, createFuncMap);

	registered = true;
      }

      return success;
    }

  }  // end namespace FormfacManageEnv

} // namespace FF
