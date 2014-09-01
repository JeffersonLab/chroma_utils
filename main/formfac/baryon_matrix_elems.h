// -*- C++ -*-
// $Id: baryon_matrix_elems.h,v 2.0 2008/12/05 04:43:47 edwards Exp $
/*! \file
 * \brief Baryon matrix elements
 */

#ifndef __baryon_matrix_elems_h__
#define __baryon_matrix_elems_h__

#include "matrix_elems.h"

namespace FF
{

  //----------------------------------------------------------------------------------
  //! Base class for meson-meson Minkowski-space matrix elements
  class BaryonMatrixElement : public MatrixElement<SpinVector,SpinVector>
  {
  public:  
    // Typedefs to save typing
    typedef SpinVector  T1;
    typedef SpinVector  T2;

    //! Virtual destructor
    virtual ~BaryonMatrixElement() {}
  };


  //! Reader
  void read(XMLReader& xml, const std::string& path, MatrixElementParams& param);

  //! Writer
  void write(XMLWriter& xml, const std::string& path, const MatrixElementParams& param);

} // namespace FF

#endif
