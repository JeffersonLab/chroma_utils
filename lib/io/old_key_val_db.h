// -*- C++ -*-
// $Id: old_key_val_db.h,v 1.2 2009/03/09 02:41:17 edwards Exp $
/*! \file
 * \brief Key and values for DB for use with ffdb-lite
 */

#ifndef __old_key_val_db_h__
#define __old_key_val_db_h__

#include "ensem/ensem.h"

#include "io/adat_io.h"
#include "DBKey.h"
#include "DBData.h"

namespace ADATIO
{
  using namespace ENSEM;
//  using namespace FFDB;

  //---------------------------------------------------------------------
  //! Serializable key harness
  /*! \ingroup ferm */
  template<typename K>
  class OldSerialDBKey : public FFDB::DBKey
  {
  public:
    //! Default constructor
    OldSerialDBKey() {} 

    //! Constructor from data
    OldSerialDBKey(const K& k) : key_(k) {}

    //! Setter
    K& key() {return key_;}

    //! Getter
    const K& key() const {return key_;}

    // Part of Serializable
    const unsigned short serialID (void) const {return 456;}

    void writeObject (std::string& output) throw (FFDB::SerializeException) {
      BinaryBufferWriter bin;
      write(bin, key());
      output = bin.str();
    }

    void readObject (const std::string& input) throw (FFDB::SerializeException) {
      BinaryBufferReader bin(input);
      read(bin, key());
    }

    // Part of DBKey
    int hasHashFunc (void) const {return 0;}
    int hasCompareFunc (void) const {return 0;}

    /**
     * Empty hash and compare functions. We are using default functions.
     */
    static unsigned int hash (Db *db, const void* bytes, unsigned int len) {return 0;}
    static int compare (Db *db, const Dbt* k1, const Dbt* k2) {return 0;}
   
  private:
    K  key_;
  };


  //---------------------------------------------------------------------
  //! Serializable value harness
  /*! \ingroup ferm */
  template<typename D>
  class OldSerialDBData : public FFDB::DBData
  {
  public:
    //! Default constructor
    OldSerialDBData() {} 

    //! Constructor from data
    OldSerialDBData(const D& d) : data_(d) {}

    //! Setter
    D& data() {return data_;}

    //! Getter
    const D& data() const {return data_;}

    // Part of Serializable
    const unsigned short serialID (void) const {return 123;}

    void writeObject (std::string& output) throw (FFDB::SerializeException) {
      BinaryBufferWriter bin;
      write(bin, data());
      ADATUtil::n_uint32_t chk(bin.getChecksum());
      write(bin, chk);
      output = bin.str();
    }

    void readObject (const std::string& input) throw (FFDB::SerializeException) {
      BinaryBufferReader bin(input);
      read(bin, data());
      ADATUtil::n_uint32_t chkd(bin.getChecksum());
      ADATUtil::n_uint32_t chk;
      read(bin, chk);
      if (! bin.eof())
      {
	if (chk != chkd)
	{
	  std::cerr << __func__ << ": data checksum mismatch\n";
	  exit(1);
	}
      }
    }

  private:
    D  data_;
  };

} // namespace ColorVec

#endif
