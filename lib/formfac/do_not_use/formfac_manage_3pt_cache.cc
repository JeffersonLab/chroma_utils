// $Id: formfac_manage_3pt_cache.cc,v 2.0 2008/12/05 04:43:35 edwards Exp $
//
// Manage building-blocks


#include "formfac/formfac_manage_3pt_cache.h"
#include "io/adat_byteorder.h"
#include <sys/time.h>
#include <string.h>

extern "C"
{
#include <lime.h>
}

#undef FF_DEBUG

namespace FF
{

#if 0
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
#endif


  //-----------------------------------------------------------------------------------
  //! Read a lime record into xml
  bool read_cache_entry(XMLReader& xml, LimeReader *reader)
  {
    n_uint64_t nbytes = limeReaderBytes(reader); // size of record
    char* xml_buf = new char[nbytes+1]; // buffer array

    n_uint64_t read_bytes = nbytes;
    int status = limeReaderReadData((void *)xml_buf, &read_bytes, reader); // read record
    if( status < 0 )
      if( status != LIME_EOR ) 
      { 
	std::cerr << __func__ << ": LIME Read Error Occurred: status= " << 
	  status << nbytes << "bytes wanted" << read_bytes 
		  << "read" << std::endl;
	delete[] xml_buf;
	return true;
      }

    // Use char buffer to initialize XMLReader
    try
    {
      std::istringstream ss;
      xml_buf[nbytes] = '\0';
      std::string foo = xml_buf;
      ss.str(foo);
      xml.open(ss);
    }
    catch (const std::string& e)
    {
      std::cerr << __func__ << ": error reading xml header: " << e << std::endl;
      return true;
    }
    
    delete[] xml_buf;

    return false;
  }


  //! Read a lime record into xml
  bool read_cache_entry(EnsemVectorComplexF& d, LimeReader *reader)
  {
    n_uint64_t nbytes = limeReaderBytes(reader); // size of record
    char* data_buf = new char[nbytes+1]; // buffer array

    n_uint64_t read_bytes = nbytes;
    int status = limeReaderReadData((void *)data_buf, &read_bytes, reader); // read record
    if( status < 0 )
      if( status != LIME_EOR ) 
      { 
	std::cerr << __func__ << ": LIME Read Error Occurred: status= " << 
	  status << nbytes << "bytes wanted" << read_bytes 
		  << "read" << std::endl;
	return true;
      }

    if (read_bytes != sizeof(ComplexF)*d.numElem()*d.size())
    {
      std::cerr << __func__ << ": lime record incorrect size to hold EnsemVectorComplexF\n";
      return true;
    }

    // Possibly byte swap
    size_t  size   = sizeof(WordType<ComplexF>::Type_t);
    int    nmemb   = sizeof(ComplexF) / sizeof(WordType<ComplexF>::Type_t);
    if (! ADATUtil::big_endian())
      ADATUtil::byte_swap(data_buf, size, read_bytes/size);

    // Snarf into each element
    d.checkSize(__func__);
    char *data_tmp = data_buf;

    for(int i=0; i < d.size(); ++i) 
    {
      d.elem(i).checkSize(__func__);
      memcpy((char*)(d.elem(i).getF()), data_tmp, d.numElem()*size*nmemb);
      data_tmp += d.numElem()*size*nmemb;
    }

    return false;
  }


  //! Write out a XML record
  void write_cache_entry(LimeWriter *writer, int MB_flag, int ME_flag,
			 XMLBufferWriter& xml, char lime_type[])
  {
    n_uint64_t bytes = strlen(xml.str().c_str())+1;

    /* Write record header */
    LimeRecordHeader* h = limeCreateHeader(MB_flag, ME_flag, lime_type, bytes);
    int status;
    if( (status = limeWriteRecordHeader( h, writer )) < 0 ) 
    { 
      std::cerr << __func__ << "L LIME write header error= " << status << std::endl;
      exit(1);
    }

    limeDestroyHeader(h);

    char *data_buf = new char[bytes];
    memcpy(data_buf, xml.str().c_str(), bytes);  // yuk, could avoid with proper const in writer below
    if( (status = limeWriteRecordData(data_buf, &bytes, writer)) != LIME_SUCCESS ) 
    { 
      std::cerr << __func__ << ": LIME write error= " << status << std::endl;
      exit(1);
    }

    delete[] data_buf;
  }


  //! Write out an Ensemble record
  void write_cache_entry(LimeWriter *writer, int MB_flag, int ME_flag,
			 const EnsemVectorComplexF& d, char lime_type[])
  {
    n_uint64_t bytes = sizeof(ComplexF)*d.numElem()*d.size();

    // Write record header
    LimeRecordHeader* h = limeCreateHeader(MB_flag, ME_flag, lime_type, bytes);
    int status;
    if( (status = limeWriteRecordHeader( h, writer )) < 0 ) 
    { 
      std::cerr << __func__ << "L LIME write header error= " << status << std::endl;
      exit(1);
    }

    limeDestroyHeader(h);

    // Snarf into each element
    char *data_buf = new char[bytes];
    char *data_tmp = data_buf;
    size_t  size   = sizeof(WordType<ComplexF>::Type_t);
    int    nmemb   = sizeof(ComplexF) / sizeof(WordType<ComplexF>::Type_t);

    for(int i=0; i < d.size(); ++i) 
    {
      memcpy(data_tmp, (char*)(d.elem(i).getF()), d.numElem()*size*nmemb);
      data_tmp += d.numElem()*size*nmemb;
    }

    // Possibly byte swap
    if (! ADATUtil::big_endian())
      ADATUtil::byte_swap(data_buf, size, bytes/size);

    if( (status = limeWriteRecordData(data_buf, &bytes, writer)) != LIME_SUCCESS ) 
    { 
      std::cerr << __func__ << ": LIME write error= " << status << std::endl;
      exit(1);
    }
  }


  //-----------------------------------------------------------------------------------
  //! Constructor
  Manage3PtFuncCache::Manage3PtFuncCache(const std::string& cache_file_,
					 int max_map_mb_) :
    cache_file(cache_file_), max_map_mb(max_map_mb_)
  {
#if 0
    // Resize the ordered_map
    std::cout << __func__ 
	      << ": bucket_count= " << threept.bucket_count()
	      << ": max_bucket_count= " << threept.max_bucket_count()
	      << ": max_load_factor= " << threept.max_load_factor() 
	      << ": max_size= " << threept.max_size() 
	      << ": size= " << threept.size() 
	      << std::endl;

    std::cout << "rehashing" << std::endl;

///    threept.rehash(20000);

    std::cout << __func__ 
	      << ": new bucket_count= " << threept.bucket_count()
	      << ": max_bucket_count= " << threept.max_bucket_count()
	      << ": max_load_factor= " << threept.max_load_factor() 
	      << ": max_size= " << threept.max_size() 
	      << ": size= " << threept.size() 
	      << std::endl;
#endif


    // Read cache file
    cache_file_exists = false;
    if (read_cache())
    {
      std::cerr << __func__ << ": fatal error in reading cache file" << std::endl;
      exit(1);
    }
  }


  //! Destructor
  Manage3PtFuncCache::~Manage3PtFuncCache()
  {
    if (cache_file_exists)
      sync_cache();
    else
      write_cache();
  }



  //! Initialize ensemble info
  void Manage3PtFuncCache::initEnsemble(const EnsembleInfo& info)
  {
    // Lots of sanity checks
    if (info.nbin <= 0)
    {
      std::cerr << __func__ << ": invalid ensemble size" << std::endl;
      exit(1);
    }

    if (info.lattice.latt_size.size() == 0)
    {
      std::cerr << __func__ << ": invalid lattice_size" << std::endl;
      exit(1);
    }

    for(int i=0; i < info.lattice.latt_size.size(); ++i)
    {
      if (info.lattice.latt_size[i] <= 0)
      {
	std::cerr << __func__ << ": invalid lattice_size extent = " << info.lattice.latt_size[i] << std::endl;
	exit(1);
      }
    }

    if (info.lattice.decay_dir < 0 || info.lattice.decay_dir > info.lattice.latt_size.size())
    {
      std::cerr << __func__ << ": invalid decay direction" << std::endl;
      exit(1);
    }

    // Finally set the info
    ensemble_info = info;
  }


  //! Read cache file
  bool Manage3PtFuncCache::read_cache()
  {
    std::cout << __func__ << std::endl;

    StopWatch swatch;
    swatch.reset();
    swatch.start();

    // Begin by opening the file
    FILE *fp;
    if((fp = fopen(cache_file.c_str(), "rb")) == NULL)
    {
      swatch.stop();
      std::cerr << "Did not find cache file " << cache_file << std::endl;
      cache_file_exists = false;
      return false;
    }

    std::cout << __func__ << ": found cache file " << cache_file << std::endl;
    cache_file_exists = true;

    // Now allocate a lime reader to read in the data
    LimeReader *reader = limeCreateReader(fp);
    if(reader == (LimeReader *) NULL)
    {
      swatch.stop();
      std::cerr << __func__ << ": unable to open LimeReader " << std::endl;
      return true;
    }

    // Move to first record
    int status;
    if ( (status = limeReaderNextRecord(reader)) != LIME_SUCCESS )
    { 
      swatch.stop();
      std::cerr << __func__ << ": limeReaderNextRecord returned status= " << status << std::endl;
      return true;
    }

    // Define a bucket to contain the XML data. 
    {
      XMLReader file_xml;
      if (read_cache_entry(file_xml, reader))
      {
	swatch.stop();
	return true;
      }

      // Extract ensemble ensemble_info
      try
      {
	EnsembleInfo cache_xml;
	read(file_xml, "/Ensem", cache_xml);

	// If the ensemble is already initialized, check that it agrees
	// with cache file
	if (isInitP())
	{
	  if (cache_xml.nbin != this->size())
	  {
	    swatch.stop();
	    std::cerr << __func__ << ": cache file has inconsistent nbins= " 
		      << cache_xml.nbin << std::endl;
	    return true;
	  }
	}
	else
	{
	  this->initEnsemble(cache_xml);
	}
      }
      catch (const std::string& e)
      {
	swatch.stop();
	std::cerr << __func__ << ": error reading xml header: " << e << std::endl;
	return true;
      }
    }


    //
    // Loop over record pairs
    //
    std::cout << __func__ << ": read records" << std::endl;

    while ((status = limeReaderNextRecord(reader)) != LIME_EOF)
    {
//    std::cout << __func__ << ": found a record" << std::endl;

      if ( status != LIME_SUCCESS ) 
      { 
	swatch.stop();
	std::cerr << __func__ << ": limeReaderNextRecord returned status= " << status << std::endl;
	return true;
      }

      n_uint64_t nbytes = limeReaderBytes(reader);
      char *lime_type   = limeReaderType(reader);
      size_t bytes_pad  = limeReaderPadBytes(reader);
      int MB_flag       = limeReaderMBFlag(reader);
      int ME_flag       = limeReaderMEFlag(reader);
    
      // Allocate space for the xml record
      XMLReader record_xml;
      if (read_cache_entry(record_xml, reader))
      {
	swatch.stop();
	return true;
      }

      // Extract three-pt info
      ThreePtArg  arg;
      try
      {
	read(record_xml, "/ThreePtArg", arg);
      }
      catch (const std::string& e)
      {
	swatch.stop();
	std::cerr << __func__ << ": error reading threept xml header: " << e << std::endl;
	return true;
      }

      // Move to next record
      if ( (status = limeReaderNextRecord(reader)) != LIME_SUCCESS )
      { 
	swatch.stop();
	std::cerr << __func__ << ": limeReaderNextRecord returned status= " 
		  << status << std::endl;
	return true;
      }

      // Snarf into ensem 
      CacheEntry_t ent;
      ent.data.resize(size());
      ent.data.resizeObs(timeLen());
      if (read_cache_entry(ent.data, reader))
      {
	swatch.stop();
	return true;
      }

      // Insert
      ent.status = CACHE_ENSEM_FILE;
      ent.count  = 1;
      threept.insert(std::make_pair(arg,ent));

    } // end while

    limeDestroyReader(reader);
    fclose(fp);

    swatch.stop();
    std::cout << __func__ << ": finished reading cache file " << cache_file 
	      << "   : total time = "
	      << swatch.getTimeInSeconds() 
	      << " secs" << std::endl;

    return false;
  }


  //! Synchronize memory with cache file
  void Manage3PtFuncCache::sync_cache()
  {
    std::cout << "Updating cache file " << cache_file << std::endl;

    StopWatch swatch;
    swatch.reset();
    swatch.start();

    if (! isInitP())
    {
      std::cerr << __func__ << ": ensemble is not initialized, not syncing" << std::endl;
      exit(1);
    }

    // Begin by opening the file
    FILE *fp;
    if((fp = fopen(cache_file.c_str(), "ab")) == NULL)
    {
      std::cerr << __func__ << ": failed to append to cache file " << cache_file << std::endl;
      exit(1);
    }

    // Now allocate a lime writer to write out the data
    LimeWriter* writer = limeCreateWriter(fp);
    if( writer == (LimeWriter *)NULL )  
    {
      std::cerr << __func__ << ": unable to open LimeWriter " << std::endl;
      exit(1);
    }

    // Loop over all marked cache entries
    // Write ones in memory but not on disk
    for(ThreePtType_t::const_iterator mm = threept.begin();
	mm != threept.end();
	++mm)
    {
      if (mm->second.status == CACHE_ENSEM_MEMORY)
      {
	XMLBufferWriter xml;
	write(xml, "ThreePtArg", mm->first);

	write_cache_entry(writer, 0, 0, xml, "XML Header");
	write_cache_entry(writer, 0, 0, mm->second.data, "Ensemble");
      }
    }

    limeDestroyWriter(writer);
    fclose(fp);

    swatch.stop();
    std::cout << __func__ << ": finished updating cache file " << cache_file 
	      << "   : total time = "
	      << swatch.getTimeInSeconds() 
	      << " secs" << std::endl;
  }



  //! Overwrite cache file
  void Manage3PtFuncCache::write_cache()
  {
    std::cout << "Writing cache file " << cache_file << std::endl;

    StopWatch swatch;
    swatch.reset();
    swatch.start();

    if (! isInitP())
    {
      std::cerr << __func__ << ": ensemble is not initialized, not writing" << std::endl;
      exit(1);
    }

    // Begin by opening the file
    FILE *fp;
    if((fp = fopen(cache_file.c_str(), "wb")) == NULL)
    {
      std::cerr << __func__ << ": failed to open for writing cache file " 
		<< cache_file << std::endl;
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
      write(xml, "Ensem", ensemble_info);

      write_cache_entry(writer, 1, 0, xml, "Cache File Header");
    }

    // Loop over all marked cache entries
    for(ThreePtType_t::const_iterator mm = threept.begin();
	mm != threept.end();
	++mm)
    {
#if defined(USE_UNORDERED_MAP)
#if 0
      std::cout << "write: "
		<< "  bucket= " << threept.bucket(mm->first)
		<< "  bucket_size= " << threept.bucket_size(threept.bucket(mm->first))
		<< "  hash= " << ThreePtArgTraits()(mm->first) 
		<< "  "
		<< "  status= " << mm->second.status 
		<< "  count= " << mm->second.count 
		<< "  p_i= " << mm->first.pi_pf.p_i
		<< "  g= " << mm->first.g
		<< "  src= " << mm->first.src
		<< "  snk= " << mm->first.snk
		<< "  links= " << mm->first.links
		<< "\n";
#endif
#endif
     
      if (mm->second.status == CACHE_ENSEM_MEMORY || mm->second.status == CACHE_ENSEM_FILE)
      {
	XMLBufferWriter xml;
	write(xml, "ThreePtArg", mm->first);
#ifdef FF_DEBUG
	std::cout << "write: status=" << mm->second.status 
		  << "  count=" << mm->second.count 
		  << "  g=" << mm->first.g
		  << "  src=" << mm->first.src
		  << "  snk=" << mm->first.snk
		  << "  links=" << mm->first.links
		  << "\n";
#endif

	write_cache_entry(writer, 0, 0, xml, "XML Header");
	write_cache_entry(writer, 0, 0, mm->second.data, "Ensemble");
      }
    }

    limeDestroyWriter(writer);
    fclose(fp);

    swatch.stop();
    std::cout << __func__ << ": finished writing cache file " << cache_file 
	      << "   : total time = "
	      << swatch.getTimeInSeconds() 
	      << " secs" << std::endl;
  }


  //! Read files
  EnsemVectorComplexF Manage3PtFuncCache::operator[](const ThreePtArg& arg)
  {
    // Fetch if not found
    if (threept.find(arg) == threept.end())
      do_read3pt(arg);

    // Mark as actually used
    ThreePtType_t::iterator it = threept.find(arg);
    if (it != threept.end())
    {
      if (it->second.status != CACHE_ENSEM_FILE)
	it->second.status = CACHE_ENSEM_MEMORY;   // mark as only in memory and used

      it->second.count++;
    }
    else
    {
      std::cerr << "operator[] internal error: entry not found" << std::endl;
      exit(1);
    }

#ifdef FF_DEBUG
    std::cout << "operator[]: status=" << it->second.status 
	      << "  count=" << it->second.count 
	      << "  g=" << it->first.g
	      << "  src=" << it->first.src
	      << "  snk=" << it->first.snk
	      << "  links=" << it->first.links
	      << "\n";
#endif

    return it->second.data;
  }
  

  //! Return map size
  int Manage3PtFuncCache::map_size_megabytes() const
  {
    // Do this the dumb way and simply count all the entries - could be slow...
    // Otherwise, I have no way of knowing if insertions are duplicates
    size_t cnt = 0;
    for(ThreePtType_t::const_iterator mm = threept.begin();
	mm != threept.end();
	++mm)
    {
      cnt++;
    }

    size_t sz = cnt*size()*timeLen()*sizeof(WordType<ComplexF>::Type_t);
    int szz = sz >> 20;
#ifdef FF_DEBUG
    std::cout << "map_size_megabytes = " << szz 
	      << "  num_map= " << cnt << "  sz=" << sz << "\n";
#endif
    return szz;
  }


  // Insert a key,value
  void Manage3PtFuncCache::insert(const ThreePtArg& arg, const EnsemVectorComplexF& v)
  {
    if (! isInitP())
    {
      std::cerr << __func__ << ": ensemble not initialized" << std::endl;
      exit(1);
    }

    if (v.numElem() != timeLen())
    {
      std::cerr << __func__ << ": inconsistent length of correlators" << std::endl;
      exit(1);
    }

    if (v.size() != size())
    {
      std::cerr << __func__ << ": inconsistent length of bins" << std::endl;
      exit(1);
    }

    // If map is filled already, then free some if this feature is turned on.
    if (max_map_mb > 0)
    {
      if (map_size_megabytes() > max_map_mb)
      {
	for(ThreePtType_t::const_iterator mm = threept.begin();
	    mm != threept.end();
	    ++mm)
	{
	  if (mm->second.status == CACHE_NOT_USED)
	  {
#ifdef FF_DEBUG
	    std::cout << "operator[]: erase; status=" << mm->second.status 
		      << "  count=" << mm->second.count 
		      << "  g=" << mm->first.g
		      << "  src=" << mm->first.src
		      << "  snk=" << mm->first.snk
		      << "  links=" << mm->first.links
		      << "\n";
#endif
	    
	    threept.erase(mm->first);
	    break;
	  }
	}
      }
    }

    // Insert
    CacheEntry_t ent;
    ent.data   = v;
    ent.status = CACHE_NOT_USED;
    ent.count  = 0;

#ifdef FF_DEBUG
    std::cout << "insert: status=" << ent.status 
	      << "  count=" << ent.count 
	      << "  g=" << arg.g
	      << "  src=" << arg.src
	      << "  snk=" << arg.snk
	      << "  links=" << arg.links
	      << "\n";
#endif

    threept.insert(std::make_pair(arg,ent));
  }


  // Does key exist?
  bool Manage3PtFuncCache::exist(const ThreePtArg& arg)
  {
    return (threept.find(arg) == threept.end()) ? false : true;
  }


  // Erase an entry
  void Manage3PtFuncCache::erase(const ThreePtArg& arg)
  {
    ThreePtType_t::iterator it = threept.find(arg);
    if (it != threept.end())
    {
      threept.erase(it);
    }
  }


  // Erase unused entries
  int Manage3PtFuncCache::eraseUnused()
  {
    int cnt = 0;

    // To avoid messing the map class, build up a list of things to erase
    // and then run through the list
    std::list<ThreePtArg> arg_list;

    // Build the list
    for(ThreePtType_t::const_iterator mm = threept.begin();
	mm != threept.end();
	++mm)
    {
      if (mm->second.status == CACHE_NOT_USED)
      {
#ifdef FF_DEBUG
	std::cout << __func__ << ": status=" << mm->second.status 
		  << "  count=" << mm->second.count 
		  << "  g=" << mm->first.g
		  << "  src=" << mm->first.src
		  << "  snk=" << mm->first.snk
		  << "  links=" << mm->first.links
		  << "\n";
#endif

	arg_list.push_back(mm->first);
      }
    }

    // Now erase the list
    for(std::list<ThreePtArg>::const_iterator mm = arg_list.begin();
	mm != arg_list.end();
	++mm)
    {
      threept.erase(*mm);
      cnt++;
    }
    return cnt;
  }


  // Mark an entry
  bool Manage3PtFuncCache::mark(const ThreePtArg& arg)
  {
    ThreePtType_t::iterator mm = threept.find(arg);
    if (mm != threept.end())
    {
#ifdef FF_DEBUG
      std::cout << __func__ << ": status=" << mm->second.status 
		<< "  count=" << mm->second.count 
		<< "  g=" << mm->first.g
		<< "  src=" << mm->first.src
		<< "  snk=" << mm->first.snk
		<< "  links=" << mm->first.links
		<< "\n";
#endif
      if (mm->second.status == CACHE_NOT_USED)
      {
	mm->second.status = CACHE_ENSEM_MEMORY;
	mm->second.count++;
      }
      return true;
    }
    else
    {
      return false;
    }
  }

} // namespace FF
