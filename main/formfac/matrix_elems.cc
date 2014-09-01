// $Id: matrix_elems.cc,v 2.0 2008/12/05 04:43:48 edwards Exp $
/*! \file
 * \brief Matrix elements
 */

#include "matrix_elems.h"

namespace FF
{
  //---------------------------------------------------------------------------
  // Read parameters
  void read(XMLReader& xml, const std::string& path, MatrixElementParams& param)
  {
    MatrixElementParams tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const std::string& path, const MatrixElementParams& param)
  {
    param.writeXML(xml, path);
  }


  //! Initialize
  MatrixElementParams::MatrixElementParams()
  {
  }


  //! Read parameters
  MatrixElementParams::MatrixElementParams(XMLReader& xml, const std::string& path)
  {
    XMLReader paramtop(xml, path);

//    read(paramtop, "MatrixElementType", mat_elem_type);
    read(paramtop, "LatticeParam", lattice);
  }

  //! Write parameters
  void MatrixElementParams::writeXML(XMLWriter& xml, const std::string& path) const
  {
    push(xml, path);

//      write(xml, "MatrixElementType",  matrix_elem_type);
    write(xml, "LatticeParam", lattice);

    pop(xml);
  }



} // namespace FF
