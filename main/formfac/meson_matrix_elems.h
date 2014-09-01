// -*- C++ -*-
// $Id: meson_matrix_elems.h,v 2.0 2008/12/05 04:43:48 edwards Exp $
/*! \file
 * \brief Meson matrix elements
 */

#ifndef __meson_matrix_elems_h__
#define __meson_matrix_elems_h__

#include "matrix_elems.h"

namespace FF
{

  //----------------------------------------------------------------------------------
  //! Base class for meson-meson Minkowski-space matrix elements
  class MesonMatrixElement : public MatrixElement<Complex,Complex>
  {
  public:  
    // Typedefs to save typing
    typedef Complex  T1;
    typedef Complex  T2;

    //! Virtual destructor
    virtual ~MesonMatrixElement() {}
  };


  //----------------------------------------------------------------------------------
  //! Pion-V-pion matrix element
  /*!
   * Matrix element is  \f$<\pi_0|V_\mu|\pi_0>\f$
   */
  class Pi0VecCurPi0MatElem : public MesonMatrixElement
  {
  public:
    // Typedefs to save typing
    typedef Complex  T1;
    typedef Complex  T2;

    //! Full constructor
    Pi0VecCurPi0MatElem(const MatrixElementParams& p) : params(p) {}

    //! Destructor
    ~Pi0VecCurPi0MatElem() {}

    //! The source state of this matrix element
    std::string sourceName() const {return "pi0";}

    //! The sink state of this matrix element
    std::string sinkName() const {return "pi0";}

    //! The number of source polarizations
    int numSourcePolar() const {return 0;}

    //! The number of sink polarizations
    int numSinkPolar() const {return 0;}

    //! The number of directions
    int numDir() const {return Nd;}

    //! The number of form-factors
    int numFF() const {return 1;}

    //! Return names for the form-factors
    std::vector<std::string> getFFNames() const;

    //! The matrix element value
    MatElemRes_t<T1,T2> operator()(const std::string& FF_name,
				   const Real& mass_f, const ArrayInt& p_f, 
				   const Real& mass_i, const ArrayInt& p_i, 
				   const Array<int>& lorentz, int r_f, int r_i) const;

  private:
    MatrixElementParams  params;   /*!< matrix element params */
  };


  //! Reader
  void read(XMLReader& xml, const std::string& path, MatrixElementParams& param);

  //! Writer
  void write(XMLWriter& xml, const std::string& path, const MatrixElementParams& param);

} // namespace FF

#endif
