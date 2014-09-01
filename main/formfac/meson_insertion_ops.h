// -*- C++ -*-
// $Id: meson_insertion_ops.h,v 2.0 2008/12/05 04:43:48 edwards Exp $
/*! \file
 * \brief Meson insertion operators
 */

#ifndef __meson_insertion_ops_h__
#define __meson_insertion_ops_h__

#include "insertion_ops.h"
#include "meson_matrix_elems.h"
#include "current_ops.h"

namespace FF
{

  //----------------------------------------------------------------------------------
  //! Base class for meson-meson Minkowski-space matrix elements
  class MesonInsertionOperator : public InsertionOperator<EnsemComplex,EnsemComplex>
  {
  public:  
    //! Virtual destructor
    virtual ~MesonInsertionOperator() {}

    //! The matrix element for this insertion
    /*! Downcast by covariant return rule */
    virtual const EnsemMesonMatrixElement& matrixElement() const = 0;

    //! The current for this insertion
    /*! Downcast by covariant return rule */
    virtual const CurrentOperator& currentOperator() const = 0;
  };


  //----------------------------------------------------------------------------------
  //! Pion-V-pion insertion
  /*!
   * Matrix element is  \f$<\pi_0|V_\mu|\pi_0>\f$
   */
  class Pi0VecCurPi0InsertOp : public MesonInsertionOperator
  {
  public:
    //! Full constructor
    Pi0VecCurPi0InsertOp(const InsertionOperatorParams& p);

    //! Destructor
    ~Pi0VecCurPi0InsertOp() {}

    //! The matrix element for this insertion
    /*! Downcast by covariant return rule */
    const EnsemMesonMatrixElement& matrixElement() const {return *mat_elem;}

    //! The current for this insertion
    /*! Downcast by covariant return rule */
    const CurrentOperator& currentOperator() const {return *current;}

  private:
    InsertionOperatorParams          params;   /*!< matrix element params */
    Handle<EnsemMesonMatrixElement>  mat_elem;
    Handle<CurrentOperator>          current;
  };



} // namespace FF

#endif
