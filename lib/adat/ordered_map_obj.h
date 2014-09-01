// -*- C++ -*-
// $Id: ordered_map_obj.h,v 1.1 2009/05/19 15:50:44 dudek Exp $
/*! \file
 * \brief Wrapper over ordered maps
 */

#ifndef __ordered_map_obj_h__
#define __ordered_map_obj_h__

#include <map>
#include <vector>

#warning "Include ordered_map"


namespace ADAT
{

  //----------------------------------------------------------------------------
  //! A wrapper over maps
  template<typename K, typename V, typename L>
  class OrderedMapObject
  {
  public:
    //! Map type convenience
    typedef std::map<K,V,L> MapType_t;

    //! Default constructor
    OrderedMapObject() {}

    //! Destructor
    ~OrderedMapObject() {}

    //! Exists?
    bool exist(const K& key) const {
      return (src_map.find(key) == src_map.end()) ? false : true;
    }
			
    //! Clear the object
    void clear() {src_map.clear();}

    //! Erase a key-value
    void erase(const K& key) {
      ADATIO::BinaryBufferWriter bin;
      write(bin, key);

      typename MapType_t::const_iterator iter = src_map.find(bin.str());
      if (iter != src_map.end())
      {
	src_map.erase(iter);
      }
    }

    //! Insert
    void insert(const K& key, const V& val) {
      //Added this to handle overwriting existing keys
      if( exist(key) )
      {
	src_map.erase(src_map.find(key));
      }

      src_map.insert(std::make_pair(key,val));
    }
			
    //! Accessor
    const V& operator[](const K& key) const {
      if (! exist(key) )
      {
	std::cerr << "MapObject: key not found" << std::endl;
	exit(1);
      }

      return src_map.find(key)->second;
    }
			
    //! The number of elements
    typename MapType_t::size_type size() const {return src_map.size();}

    //! Dump keys
    std::vector<K> keys() const {
      std::vector<K> keys;
      typename MapType_t::const_iterator iter;
      for(iter  = src_map.begin();
	  iter != src_map.end();
	  ++iter)
      {
	keys.push_back(iter->first);
      }
      return keys;
    }

    //! Usual begin iterator
    typename MapType_t::const_iterator begin() const {return src_map.begin();}

    //! Usual end iterator
    typename MapType_t::const_iterator end() const {return src_map.end();}

  private:
    //! Map of objects
    mutable MapType_t src_map;
  };

} // namespace ColorVec

#endif
