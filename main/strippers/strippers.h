// -*- C++ -*-
// $Id: strippers.h,v 2.0 2008/12/05 04:44:05 edwards Exp $

/*! \file
 * \brief Some shortcuts
 */


#ifndef __strippers_h__
#define __strippers_h__

#include "xml_array.h"
#include "io/adat_xmlio.h"
#include "io/adat_io.h"

namespace Strippers 
{
  // namespace composition
  using XMLArray::Array;
  using namespace ADATXML;
  using namespace ADATIO;


  // Types snarfed from ADAT. These should be merged somehow
  typedef float  Real;
  typedef double Double;
  typedef int    Integer;
  typedef bool   Boolean;

  struct Complex
  {
    Real  re;
    Real  im;
  };

  struct ComplexD
  {
    Double  re;
    Double  im;
  };


  // These readers are here in the Strippers namespace to facilitate
  // an Array<Complex> read. If not in this namespace, the Array read
  // will not know to where to look.
  inline
  void read(XMLReader& xml, const std::string& path, Complex& com)
  {
    try {
      XMLReader complextop(xml, path);
      read(complextop, "re", com.re);
      read(complextop, "im", com.im);
    }
    catch(const std::string& error) { 
      std::ostringstream error_stream;
      error_stream << "Error reading complex: " << error << std::endl;
      throw error_stream.str();
    }
  }


  // These readers are here in the Strippers namespace to facilitate
  // an Array<Complex> read. If not in this namespace, the Array read
  // will not know to where to look.
  inline
  void read(XMLReader& xml, const std::string& path, ComplexD& com)
  {
    try {
      XMLReader complextop(xml, path);
      read(complextop, "re", com.re);
      read(complextop, "im", com.im);
    }
    catch(const std::string& error) { 
      std::ostringstream error_stream;
      error_stream << "Error reading complex: " << error << std::endl;
      throw error_stream.str();
    }
  }


  // This needs to be in Strippers namespace, otherwise an Array<Complex>
  // read will not find it.
  inline
  void read(BinaryReader& bin, Complex& p)
  {
    read(bin, p.re);
    read(bin, p.im);
  }

}  // namespace Strippers


#endif
