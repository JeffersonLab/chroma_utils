// -*- C++ -*-
// $Id: adat_qdpio.h,v 2.0 2008/12/05 04:43:37 edwards Exp $

/*! @file
 * @brief IO support via QIO
 */

#ifndef ADAT_QDPIO_H
#define ADAT_QDPIO_H

#include "xml_array.h"
#include "io/adat_xmlio.h"
#include "ensem/ensem.h"

extern "C"
{
#include <lime.h>
}
                                                                                
namespace ADATIO
{
  using XMLArray::Array;
  using namespace ENSEM;
  using namespace ADATXML;

  //--------------------------------------------------------------------------------
  /*! @defgroup qio QIO
   *
   *
   * @{
   */

  //! File open mode
  enum ADAT_filemode_t
  {
    QDPIO_CREATE,
    QDPIO_OPEN,
    QDPIO_APPEND,
  };

  //--------------------------------------------------------------------------------
  //! QIO private file info
  struct QIOFileInfo_t
  {
    Array<int>    dims;        /*!< Lattice size */
  };


  //! QIO private record info
  struct QIORecordInfo_t
  {
    std::string    date;          /*!< Date written */
    int            recordtype;    /*!< 0 if field 1 if global 2 if hypercube subset 0 */
    std::string    datatype;      /*!< Some kind of QLA like name */
    std::string    precision;     /*!< I, F, D, or S (random no state) */
    int            colors;        /*!< number of colors */
    int            spins;         /*!< number of spins */
    int            typesize;      /*!< byte length of datum */
    int            datacount;     /*!< number of datum */
  };


  //--------------------------------------------------------------------------------
  //! QIO Type and precision strings. 
  //  Using the magic of C++ I can define the right type and precision
  //  strings i Need to pass to QIO using templates. To do this I need
  //  templated structures with static members.

  //! Catch all case
  template<typename T>
  struct QIOStringTraits 
  {
    static char* tname;
    static char* tprec;
  };

  //! Partial(?) Specialisation for EScalar Objects
  template<typename T>
  struct QIOStringTraits<EScalar<T> >
  {
    static char* tname;
    static char* tprec;
  };

  //! Partial(?) Specialisation for EScalar Objects
  template<typename T>
  struct QIOStringTraits<Array< EScalar<T> > >
  {
    static char* tname;
    static char* tprec;
  };

  //! Generic type
  template<typename T>
  char* QIOStringTraits<T>::tname = "QDP_GenericType";

  //! Scalar Type
  template<typename T> 
  char* QIOStringTraits< EScalar<T> >::tname = "Scalar";

  //! Array<ScalarType>
  template<typename T> 
  char* QIOStringTraits< Array<EScalar<T> > >::tname = "Scalar";

  //! Unknown precision string
  template<typename T>
  char*  QIOStringTraits<T>::tprec = "U"; 

  // Full specialisations deferred to the qdp_qio_strings.cc file
  template<>
  char* QIOStringTraits<float>::tprec;

  template<>
  char* QIOStringTraits<double>::tprec;

  template<>
  char* QIOStringTraits<int>::tprec;
  

  //--------------------------------------------------------------------------------
  //! QIO class
  /*!
    This is a QDP object wrapper around the QIO library.
 
    QIO is a C library independentof QDP. It is designed to read/write SCIDAC
    format data files, which means a mixture of binary data and XML
    metadata together in the same file according to a scheme called Lime.
    There is a seperate independent library for handling general Lime files.

    The data is assumed to be a record in a Lime file. The user metadata (both
    file and record) is also read.
 
    Data is assumed to be stored in the file in big-endian format and any
    necessary byte-swapping is taken care of.

    The status of the IO operations is monitored internally and can be queried.
  */

  class ADATFileReader
  {
  public:
    //! Partial constructor
    ADATFileReader();

    //! Closes the last file opened.
    ~ADATFileReader();

    //! Opens a file for reading
    /*!
      Also reads the file user metadata record. 
      File is only opened in serial mode.
      \param xml Container for the file metadata.
      \param path The name of the file.
    */
    ADATFileReader(XMLReader& xml, const std::string& path);
  
    //! Opens a file for reading
    /*!
      Also reads the file user metadata record.
      File is only opened in serial mode.
      \param xml Container for the file metadata.
      \param path The name of the file.
    */
    void open(XMLReader& xml, const std::string& path);

    //! Closes the last file opened.
    void close();

    //! Queries whether a file is open
    /*!
      \return true if a file is open; false otherwise.
    */
    bool is_open();

    //! Reads an EScalar object
    template<class T>
    void read(XMLReader& xml, EScalar<T>& s1);

    //! Reads an array of objects all in a single record
    template<class T>
    void read(XMLReader& xml, Array< EScalar<T> >& s1);

    //! Reads as a string plus BinaryBufferReader pair
    void read(std::string& xml_str, BinaryBufferReader& s1);

    //! Reads an XMLReader plus BinaryBufferReader pair
    void read(XMLReader& xml, BinaryBufferReader& s1);

    //! Query whether the end-of-file has been reached.
    /*!
      \return True if  the end-of-file has been reached; false otherwise
    */
    bool eof() const;

    //! Return the lattice size 
    /*! This will not be valid if the file has not been opened */
    Array<int> lattSize() const {return file_info.dims;}

  protected:
    //! Read one single record filled with a string
    void readXML(std::string& xml_str);

    //! Read one single record filled with XML
    void readXML(XMLReader& xml_in);

    //! Read one single record filled with binary into a char array
    /*! 
      NOTE: we will always allocate one extra character to be used by the
      driver code, such as nulls.

      \param bin        The array that will hold the data
      \param word_size  The word size in bytes of each data. Used for byte reversal
     */
    void readBinary(Array<char>& bin);

    //! Convenience function to read one logical record 
    /*! 
      NOTE: we will always allocate one extra character to be used by the
      driver code, such as nulls.
     */
    void readLogicalRecord(XMLReader& xml_in, Array<char>& bin);

    //! Get access to the lime reader
    LimeReader *get() const {return lime_reader;}

    //! Read record info
    void readRecordInfo(std::string& xml_str);

    //! Read record info
    void readRecordInfo(XMLReader& xml_in);

    //! Read binary and checksum
    void readRecordData(void* bin, int datum_size, int word_size);

    //! Move to next record
    void nextRecord();

    //! Convenience function for exiting
    void blowup(const std::string& error_string);

  private:
    bool               iop;
    QIOFileInfo_t      file_info;
    QIORecordInfo_t    record_info;
    FILE*              fp;
    LimeReader*        lime_reader;
  };


  // Convenience functions

  //! Reads an EScalar object
  /*!
    \param qsw The reader
    \param rec_xml The user record metadata.
    \param sl The data
  */
  template<class T>
  void read(ADATFileReader& qsw, XMLReader& rec_xml, EScalar<T>& s1)
  {
    qsw.read(rec_xml,s1);
  }

  //! Reads an array of EScalar object
  /*!
    \param qsw The reader
    \param rec_xml The user record metadata.
    \param sl The data
  */
  template<class T>
  void read(ADATFileReader& qsw, XMLReader& rec_xml, Array< EScalar<T> >& s1)
  {
    qsw.read(rec_xml,s1);
  }

  //! Reads a BinaryBufferReader object
  /*!
    \param qsw The reader
    \param rec_xml The user record metadata.
    \param sl The data
  */
  inline
  void read(ADATFileReader& qsw, XMLReader& rec_xml, BinaryBufferReader& s1)
  {
    qsw.read(rec_xml,s1);
  }

  //! Closes a ADATFileReader.
  void close(ADATFileReader& qsw);

  //! Queries whether a ADATFileReader is open.
  bool is_open(ADATFileReader& qsw);




  //-------------------------------------------------
  //! Reads an EScalar object
  /*!
    \param rec_xml The (user) record metadata.
    \param sl The data
  */

  template<typename T>
  void ADATFileReader::read(XMLReader& rec_xml, EScalar<T>& s1)
  {
    readRecordInfo(rec_xml);
  
    switch ( record_info.precision[0] ) 
    { 
    case 'F':
    {
      std::cout << "Single Precision Read" << std::endl;
      EScalar< typename SinglePrecType<T>::Type_t > from_disk;
      readRecordData((void *)&(from_disk.elem()),
		     sizeof(typename SinglePrecType<T>::Type_t),
		     sizeof(typename WordType< typename SinglePrecType<T>::Type_t >::Type_t));

      std::cout << "readRecordData finished" << std::endl;
      s1 = from_disk;
    }
    break;
    case 'D' :
    {
      std::cout << "Reading Double Precision" << std::endl;
      EScalar< typename DoublePrecType<T>::Type_t > from_disk;
      readRecordData((void *)&(from_disk.elem()),
		     sizeof(typename DoublePrecType<T>::Type_t),
		     sizeof(typename WordType< typename DoublePrecType<T>::Type_t >::Type_t));
      s1 = from_disk;
      std::cout << "readRecordData finished" << std::endl;
    }
    break;
    default:
    {
      std::cout << "Reading I or U Precision" << std::endl;
      readRecordData((void *)&(s1.elem()),
		     sizeof(T),
		     sizeof(typename WordType<T>::Type_t));
      std::cout << "QIO_read_finished" << std::endl;
    }
    break;
    }
  }
  

  //! Reads an array of EScalar objects
  /*!
    This implementation is only correct for scalar ILattice

    \param rec_xml The (user) record metadata.
    \param sl The data
  */
  template<typename T>
  void ADATFileReader::read(XMLReader& rec_xml, Array< EScalar<T> >& s1)
  {
    readRecordInfo(rec_xml);
  
    switch ( record_info.precision[0] ) 
    { 
    case 'F':
    {
      std::cout << "Single Precision Read" << std::endl;
      Array< EScalar< typename SinglePrecType<T>::Type_t > > from_disk(record_info.datacount);
      readRecordData((void *)from_disk.slice(),
		     from_disk.size()*sizeof(typename SinglePrecType<T>::Type_t),
		     sizeof(typename WordType< typename SinglePrecType<T>::Type_t >::Type_t));

      std::cout << "readRecordData finished" << std::endl;
      
      // Cast appropriately
      s1.resize(from_disk.size());
      for(int i=0; i < from_disk.size(); i++) { 
	s1[i] = from_disk[i];
      }
    }
    break;

    case 'D':
    {
      std::cout << "Reading Double Precision" << std::endl;
      Array< EScalar< typename DoublePrecType<T>::Type_t > > from_disk(record_info.datacount);
      readRecordData((void *)from_disk.slice(),
		     from_disk.size()*sizeof(typename DoublePrecType<T>::Type_t),
		     sizeof(typename WordType< typename DoublePrecType<T>::Type_t >::Type_t));

      std::cout << "readRecordData finished" << std::endl;
      
      // Cast appropriately
      s1.resize(from_disk.size());
      for(int i=0; i < from_disk.size(); i++) { 
	s1[i] = from_disk[i];
      }
    }
    break;

    default:
    {
      std::cout << "Reading I or U Precision" << std::endl;
      s1.resize(record_info.datacount);
      readRecordData((void *)s1.slice(),
		     s1.size()*sizeof(T),
		     sizeof(typename WordType<T>::Type_t));

      std::cout << "readRecordData finished" << std::endl;
    }
    break;
    }
  }


  // Reads a BinaryBufferReader object
  /*!
    \param rec_xml The (user) record metadata.
    \param sl The data
  */
//  void ADATFileReader::read(XMLReader& rec_xml, BinaryBufferReader& s1);


} // namespace ADATIO
#endif
