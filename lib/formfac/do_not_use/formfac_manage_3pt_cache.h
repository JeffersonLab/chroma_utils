// -*- C++ -*-
// $Id: formfac_manage_3pt_cache.h,v 2.0 2008/12/05 04:43:35 edwards Exp $
/*! \file
 * \brief Manage 3-pt funcs 
 */

#ifndef __formfac_manage_3pt_cache_h__
#define __formfac_manage_3pt_cache_h__

#include "formfac/formfac_manage_3pt.h"
#include "adat/map_traits.h"

// Unfortunately, have to include lime.h here. I only want forward decls,
// but C does not make this easy with a typedef of a struct.
extern "C"
{
#include <lime.h>
}

namespace FF
{
  using namespace ADAT;

  //----------------------------------------------------------------------------
  //! Class to hold 3-pt functions using a cache mechanism.
  /*! Still abstract - relies on  do_read3pt  to be implemented in derived class. */
  class Manage3PtFuncCache : public Manage3PtFuncMap
  {
  public:
    //! Constructor
    Manage3PtFuncCache(const std::string& cache_file_, 
		      int max_map_mb_);

    //! Destructor
    virtual ~Manage3PtFuncCache();

    //! Return a 3-pt ref. Fetch if not in memory
    virtual EnsemVectorComplexF operator[](const ThreePtArg& arg);

    //! Always insert a key,value
    virtual void insert(const ThreePtArg& k, const EnsemVectorComplexF& v);

    //! Erase an entry
    virtual void erase(const ThreePtArg& k);

    //! Erase all unused entries
    virtual int eraseUnused();

    //! Does key exist?
    virtual bool exist(const ThreePtArg& k);

    //! Mark an entry as needed
    virtual bool mark(const ThreePtArg& arg);

    //! Return map size in megabytes
    virtual int map_size_megabytes() const;

    //! Number of bins
    virtual int size() const {return ensemble_info.nbin;}

    //! Time extent
    virtual int timeLen() const {return ensemble_info.lattice.latt_size[ensemble_info.lattice.decay_dir];}

    //! Decay direction
    virtual int decayDir() const {return ensemble_info.lattice.decay_dir;}

    //! Lattice Size
    virtual const ArrayInt& lattSize() const {return ensemble_info.lattice.latt_size;}

    //! Is ensemble initialized
    virtual bool isInitP() const {return (ensemble_info.nbin == 0) ? false : true;}

  protected:
    //! Enums for cache file management
    enum CacheStatus_t
    {
      CACHE_NOT_USED,
      CACHE_ENSEM_MEMORY,
      CACHE_ENSEM_FILE
    };

    //! Memory entry
    struct CacheEntry_t
    {
      CacheEntry_t() {status=CACHE_NOT_USED;count = 0;}

      EnsemVectorComplexF   data;
      CacheStatus_t         status;
      int                   count;
    };

    //! Initialize ensemble info
    virtual void initEnsemble(const EnsembleInfo& _info);

    //! Return ensemble info
    virtual const EnsembleInfo& ensembleInfo() const {return ensemble_info;}

    //! Read building 3pt files
    virtual void do_read3pt(const ThreePtArg& arg) = 0;

    //! Read cache file
    virtual bool read_cache();

    //! Overwrite cache file
    virtual void write_cache();

    //! Synchronize memory with cache file
    virtual void sync_cache();

    //! Beginning
//    virtual const_iterator begin() const {return threept.begin();}

    //! Ending
//    virtual const_iterator end() const {return threept.end();}

  private:
    std::string           cache_file;
    bool                  cache_file_exists;
    int                   max_map_mb;
    EnsembleInfo          ensemble_info;

    //! TR1 unordered_map version (a hash table)
    typedef std::tr1::unordered_map<ThreePtArg,
				    CacheEntry_t,
				    ADAT::UnorderedMapTraits<ThreePtArg>,
				    ADAT::UnorderedMapTraits<ThreePtArg> > ThreePtType_t;

    //! Memory version of a data-base
    ThreePtType_t         threept;
  };


  //! Read a lime record into xml
  bool read_cache_entry(XMLReader& xml, LimeReader* reader);

  //! Read a lime record into xml
  bool read_cache_entry(EnsemVectorComplexF& d, LimeReader* reader);

  //! Write out a XML record
  void write_cache_entry(LimeWriter* writer, int MB_flag, int ME_flag,
			 XMLBufferWriter& xml, char lime_type[]);
  
  //! Write out an Ensemble record
  void write_cache_entry(LimeWriter* writer, int MB_flag, int ME_flag,
			 const EnsemVectorComplexF& d, char lime_type[]);


} // namespace FF

#endif
