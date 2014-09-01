// $Id: adat_qdpio.cc,v 2.0 2008/12/05 04:43:37 edwards Exp $
//
/*! @file
 * @brief IO support via QIO
 */

#include "io/adat_qdpio.h"
#include "io/adat_byteorder.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;

namespace ADATIO
{

  //! QIO record checksum
  struct QIOChecksumInfo_t
  {
    n_uint32_t    suma;        /*!< first checksum */
    n_uint32_t    sumb;        /*!< second checksum */
  };


  //! Read a private file info struct
  void read(XMLReader& xml, const std::string& path, QIOFileInfo_t& param)
  {
    XMLReader paramtop(xml, path);

    string version;
    read(paramtop, "version", version);

    int spacetime;
    read(paramtop, "spacetime", spacetime);

    if (version == "1.1")
    {
      // nop
    }
    else 
    {
      cerr << __func__ << ": unsupported version in reading QIOFileInfo_t struct: version=" 
	   << version
	   << endl;
      exit(1);
    }

    read(paramtop, "dims", param.dims);  // about the only useful thing
    if (param.dims.size() != spacetime)
    {
      cerr << __func__ << ": inconsistent dimension size in private file xml" << endl;
      exit(1);
    }

    int volfmt;
    read(paramtop, "volfmt", volfmt);  // only support singlefile
    if (volfmt != 0)
    {
      cerr << __func__ << ": only support singlefile format: volfmt=" 
	   << volfmt
	   << endl;
      exit(1);
    }
  }


  //! Read a private record info struct
  void read(XMLReader& xml, const std::string& path, QIORecordInfo_t& param)
  {
    XMLReader paramtop(xml, path);

    string version;
    read(paramtop, "version", version);

    if (version == "1.0")
    {
      int globaldata;
      read(paramtop, "globaldata", globaldata);

      if (globaldata != 1)
      {
	cerr << __func__ << ": QIORecordInfo reader: only support global data" << endl;
	exit(1);
      }
    }
    else if (version == "1.1")
    {
      int recordtype;
      read(paramtop, "recordtype", recordtype);
      
      if (recordtype != 1)
      {
	cerr << __func__ << ": QIORecordInfo reader: only support global data" << endl;
	exit(1);
      }
    }
    else 
    {
      cerr << __func__ << ": unsupported version in reading QIORecordInfo_t struct: version=" 
	   << version
	   << endl;
      exit(1);
    }
    
    read(paramtop, "date", param.date);
    read(paramtop, "datatype", param.datatype);
    read(paramtop, "precision", param.precision);
    read(paramtop, "typesize", param.typesize);
    read(paramtop, "datacount", param.datacount);
  }


  //! Private functions
  namespace 
  {
    //! Convert hex to int
    n_uint32_t convHexToInt(const string& s)
    {
      n_uint32_t v;

      std::istringstream is(s);
      is >> std::hex;
      is >> v;

      return v;
    }
  }


  //! Read a checksum info struct
  void read(XMLReader& xml, const std::string& path, QIOChecksumInfo_t& param)
  {
    XMLReader paramtop(xml, path);

    string version;
    read(paramtop, "version", version);
    if (version != "1.0")
    {
      cerr << "Unsupported version of checksum reader: version=" << version << endl;
      exit(1);
    }

    string chk;
    read(paramtop, "suma", chk); param.suma = convHexToInt(chk);
    read(paramtop, "sumb", chk); param.sumb = convHexToInt(chk);

    printf("suma= %x  sumb= %x\n", param.suma, param.sumb);
  }


  //-----------------------------------------------------------------------------
  // QDP QIO support
  ADATFileReader::ADATFileReader() {iop=false;}

  ADATFileReader::ADATFileReader(XMLReader& xml, 
				 const std::string& path)
  {open(xml,path);}

  //! Read one single record filled with XML
  void ADATFileReader::readXML(std::string& xml_str)
  {
    cout << __func__ << ": reading XML record: type = XX" 
	 << limeReaderType(lime_reader) 
	 << "XX" << endl;

    //  Allocate space for the xml record
    n_uint64_t nbytes = limeReaderBytes(lime_reader)+1; // size of record

    char* xml_buf = new char[nbytes+2]; // buffer array

    int status = limeReaderReadData((void *)xml_buf, &nbytes, lime_reader); // read record
    if ( status != 0 )
    {
      cerr <<  "LIME Read Error Occurred: status= " << status 
	   << lime_reader->bytes_total << "bytes wanted" << nbytes 
	   << "read" << endl;
      limeDestroyReader(lime_reader); 
      fclose(fp);
      exit(1);
    }

    // Use char buffer to initialize XMLReader
    xml_buf[nbytes+1] = '\0';
    xml_str = xml_buf;

    cout << "XML = XX" << xml_str << "XX" << endl;

    delete[] xml_buf;  // cleanup
  }


  //! Read one single record filled with XML
  void ADATFileReader::readXML(XMLReader& xml_in)
  {
    // String for the xml
    std::string xml_str;
    readXML(xml_str);

    // Use string to initialize XMLReader
    std::istringstream ss;
    ss.str(xml_str);
    xml_in.open(ss);

    cout << "XML = XX" << ss.str() << "XX" << endl;
  }


  //! Read one single record filled with binary into a char array
  void ADATFileReader::readBinary(Array<char>& bin)
  {
    cout << __func__ << ": reading binary record: type = XX" 
	 << limeReaderType(lime_reader) 
	 << "XX" << endl;

    n_uint64_t nbytes = limeReaderBytes(lime_reader);
    bin.resize(nbytes+1);
    int status = limeReaderReadData((void *)(bin.slice()), &nbytes, lime_reader);
    if( status != 0 ) 
    {
      cerr << "LIME read error occurred: status= " << status
	   << lime_reader->bytes_total << " wanted " << nbytes << "read"<<endl;
      exit(1);
    }
  }


  void ADATFileReader::open(XMLReader& file_xml, 
			    const std::string& path)
  {
    try
    {
      // Open fiile
      if((fp = fopen(path.c_str(), "rb")) == NULL)
      {
	cerr << "Failed to open file " << path << endl;
	exit(1);
      }

      // Now allocate a lime reader to read in the data
      lime_reader = limeCreateReader(fp);
      if (lime_reader == (LimeReader *) NULL)
      {
	cerr << "Unable to open LimeReader " << endl;
	fclose(fp);
	exit(1);
      }
   
      // Move to first record
      nextRecord();

      // Read private file xml
      cout << __func__ << ": read private file xml" << endl;
      XMLReader private_file_xml;
      readXML(private_file_xml);

      ADATIO::read(private_file_xml, "/scidacFile", file_info);

      // Move to next record
      nextRecord();

      // Read public file xml
      cout << __func__ << ": read public file xml" << endl;
      readXML(file_xml);
    }
    catch(const std::string& e) 
    {
      cerr << __func__ << ": caught exception opening file: " << e << endl;
      exit(1);
    }
    catch(std::exception& e) 
    {
      cerr << __func__ << ": caught standard library exception opening file: " << e.what() << endl;
      exit(1);
    }
    catch(...)
    {
      cerr << __func__ << ": caught generic exception opening file" << endl;
      exit(1);
    }

    iop=true;
  }


  void ADATFileReader::close()
  {
    if (is_open()) 
    {
      // close up reader
      limeDestroyReader(lime_reader);
      fclose(fp);
    }

    iop = false;
  }

  bool ADATFileReader::is_open() {return iop;}

  bool ADATFileReader::eof() const {return false;}

  ADATFileReader::~ADATFileReader() {close();}

  void ADATFileReader::nextRecord()
  {
    int status;

    cout << __func__ << ": move to next record" << endl;

    // Move to next record
    if((status = limeReaderNextRecord(lime_reader)) != LIME_SUCCESS)
    {
      cerr << "limeReaderNextRecord returned status = "<< status<<endl;
      limeDestroyReader(lime_reader); 
      fclose(fp);
      exit(1);
    }

    cout << __func__ << ": new record type = XX" 
	 << limeReaderType(lime_reader) 
	 << "XX" << endl;
  }


  void ADATFileReader::readRecordInfo(std::string& record_xml_str)
  {
    try
    {
      // Move to first record within the next logical record
      nextRecord();

      // Read private record xml
      cout << __func__ << ": read private record xml" << endl;
      XMLReader private_record_xml;
      readXML(private_record_xml);

      ADATIO::read(private_record_xml, "/scidacRecord", record_info);

      // Move to second record within logical record
      nextRecord();

      // Read user record xml
      cout << __func__ << ": read user record xml" << endl;
      readXML(record_xml_str);
    }
    catch(const std::string& e) 
    {
      cerr << __func__ << ": error reading record info: " << e << endl;
      exit(1);
    }
  }


  void ADATFileReader::readRecordInfo(XMLReader& record_xml)
  {
    try
    {
      // Read user record xml
      std::string record_xml_str;
      readRecordInfo(record_xml_str);

      std::stringstream ss(record_xml_str);
      record_xml.open(ss);
    }
    catch(const std::string& e) 
    {
      cerr << __func__ << ": error reading record info: " << e << endl;
      exit(1);
    }
  }


  void ADATFileReader::readRecordData(void* bin, int datum_size, int word_size)
  {
    try
    {
      // Move to binary payload within logical record
      nextRecord();

      // Read record binary payload
      n_uint64_t nbytes = limeReaderBytes(lime_reader);
      if (nbytes != datum_size*word_size)
      {
	cerr << __func__ << ": number of bytes not equal to datum_size*word_size: "
	     << lime_reader->bytes_total << " wanted " << nbytes << endl;
	exit(1);
      }

      int status = limeReaderReadData(bin, &nbytes, lime_reader);
      if( status != 0 ) 
      {
	cerr << "LIME read error occurred: status= " << status
	     << lime_reader->bytes_total << " wanted " << nbytes << "read"<<endl;
	exit(1);
      }

      // Compute checksum. 
      // For global data, this is a straightforward crc32 on the entire buffer
      n_uint32_t zro = 0;
      n_uint32_t chk = ADATUtil::crc32(zro, (unsigned char *)bin, (size_t)nbytes);

      // Move to the checksum record
      nextRecord();

      // Read the checksum
      QIOChecksumInfo_t  checksum_info;
      XMLReader checksum_xml;
      readXML(checksum_xml);
      ADATIO::read(checksum_xml, "/scidacChecksum", checksum_info);

      // Check the checksums. For global data, both parts of the checksum
      // in the record are the same
      if (chk != checksum_info.suma)
      {
	cerr  << __func__ << ": fatal error - checksum mismatch" << endl;
	exit(1);
      }
      else
      {
	cout << __func__ << ": checksum okay!" << endl;
      }
    }
    catch(const std::string& e) 
    {
      cerr << __func__ << ": error reading binary payload or checksum: " << e << endl;
      exit(1);
    }
  }


  void ADATFileReader::blowup(const std::string& error_string)
  {
    cerr << error_string << endl;
    limeDestroyReader(lime_reader); 
    fclose(fp);
    exit(1);
  }


  // Convenience function to read one logical record 
  /* 
    NOTE: we will always allocate one extra character to be used by the
    driver code, such as nulls.
  */
  void ADATFileReader::readLogicalRecord(XMLReader& xml_in, Array<char>& bin)
  {
    cout << __func__ << endl;

    // Move to next record
    for(int rec=0; rec < 2; ++rec)
    {
      cout << __func__ << ": skipping record: type = XX" 
	   << limeReaderType(lime_reader) 
	   << "XX" << endl;

      int status;
      if((status = limeReaderNextRecord(lime_reader)) != LIME_SUCCESS)
      {
	cerr << "limeReaderNextRecord returned status = "<< status<<endl;
	exit(1);
      }
    }

    cout << __func__ << ": reading record: type = XX" 
	 << limeReaderType(lime_reader) 
	 << "XX" << endl;

    int status;

    //  Allocate space for the xml record
    readXML(xml_in);

    // Read the binary payload
    if ((status = limeReaderNextRecord(lime_reader)) != LIME_SUCCESS )
    {
      cerr << "limeReaderNextRecord returned status = "<< status<<endl;
      limeDestroyReader(lime_reader); 
      fclose(fp);
      exit(1);
    }

    // Read binary
    readBinary(bin);

    // Move to next record
    for(int rec=0; rec < 1; ++rec)
    {
      if((status = limeReaderNextRecord(lime_reader)) != LIME_SUCCESS)
      {
	cerr << "limeReaderNextRecord returned status = "<< status<<endl;
	limeDestroyReader(lime_reader); 
	fclose(fp);
	exit(1);
      }

      cout << __func__ << ": at end of readRecord: type = XX" 
	   << limeReaderType(lime_reader) 
	   << "XX" << endl;
    }
  }


  // Reads a BinaryBufferReader object
  /*!
    \param rec_xml The (user) record metadata.
    \param sl The data
  */
  void ADATFileReader::read(std::string& record_xml_str, BinaryBufferReader& s1)
  {
    // Read record info
    readRecordInfo(record_xml_str);

    // Read binary
    // NOTE: C++ makes it painful for me to just give the final string
    // directly into the final istringstream. If c-lime implemented
    // istreams, I'd just read the istringstream directly. Instead,
    // I must read in the string into some kind of character array,
    // then insert it in to a string, then insert the string into the
    // stringstream. Yuk.
    Array<char> from_disk(record_info.datacount+1);
    readRecordData((void *)from_disk.slice(),
		   record_info.datacount,
		   sizeof(char));

    // Move the data into a binary buffer reader
    from_disk[record_info.datacount] = '\0';
    std::string foo;
    foo.resize(record_info.datacount+1);
    foo.insert(0, from_disk.slice(), record_info.datacount+1);
    s1.open(foo);
  }


  // Reads a BinaryBufferReader object
  /*!
    \param rec_xml The (user) record metadata.
    \param sl The data
  */
  void ADATFileReader::read(XMLReader& record_xml, BinaryBufferReader& s1)
  {
    try
    {
      // Read record info
      std::string record_xml_str;
      read(record_xml_str, s1);

      std::stringstream ss(record_xml_str);
      record_xml.open(ss);
    }
    catch(const std::string& e) 
    {
      cerr << __func__ << ": error reading record: " << e << endl;
      exit(1);
    }
  }


  //! Close a ADATFileReader
  void close(ADATFileReader& qsw)
  {
    qsw.close();
  }

  //! Is a ADATFileReader open
  bool is_open(ADATFileReader& qsw)
  {
    return qsw.is_open();
  }


} // namespace ADATIO
