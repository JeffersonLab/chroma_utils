// -*- C++ -*-
// $Id: matrix_elems.h,v 2.0 2008/12/05 04:43:48 edwards Exp $
/*! \file
 * \brief Matrix elements
 */

#ifndef __matrix_elems_h__
#define __matrix_elems_h__

#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "formfac/formfac_ensemble.h"
#include <string>
#include <vector>

namespace FF
{
  // Namespace composition
  using namespace ENSEM;

  //----------------------------------------------------------------------------------
  template<typename T1, typename T2>
  struct MatElemRes_t
  {
  };


  template<>
  struct MatElemRes_t<Complex,Complex>
  {
    Complex  result;
  };


  template<>
  struct MatElemRes_t<SpinVector,Complex>
  {
    SpinVector  result;
  };


  template<>
  struct MatElemRes_t<Complex,SpinVector>
  {
    SpinVector  result;
  };


  template<>
  struct MatElemRes_t<SpinVector,SpinVector>
  {
    Complex     result;
    SpinVector  source;
    SpinVector  sink;
  };


  template<>
  struct MatElemRes_t<EnsemComplex,EnsemComplex>
  {
    EnsemComplex  result;
  };


  template<>
  struct MatElemRes_t<EnsemSpinVector,EnsemComplex>
  {
    EnsemSpinVector  result;
  };


  template<>
  struct MatElemRes_t<EnsemComplex,EnsemSpinVector>
  {
    EnsemSpinVector  result;
  };


  template<>
  struct MatElemRes_t<EnsemSpinVector,EnsemSpinVector>
  {
    EnsemComplex     result;
    EnsemSpinVector  source;
    EnsemSpinVector  sink;
  };


  //----------------------------------------------------------------------------------
  //! Base class for Minkowski-space matrix elements
  template<typename T1, typename T2>
  class MatrixElement
  {
  public:  
    //! Virtual destructor
    virtual ~MatrixElement() {}

    //! The source state of this matrix element
    virtual std::string sourceName() const = 0;

    //! The sink state of this matrix element
    virtual std::string sinkName() const = 0;

    //! The number of source polarizations
    virtual int numSourcePolar() const = 0;

    //! The number of sink polarizations
    virtual int numSinkPolar() const = 0;

    //! The number of directions in the insertion
    virtual int numDir() const = 0;

    //! The number of form-factors
    virtual int numFF() const = 0;

    //! Return names for the form-factors
    virtual std::vector<std::string> getFFNames() const = 0;

    //! The matrix element value multiplying the form-factor given by a string name
    /*! 
     *  \param FF_name  Name of the form-factor
     *  \param mass_f   Final state mass
     *  \param p_f      Minkowski-space 3-vector
     *  \param mass_i   Initial state mass
     *  \param p_i      Minkowski-space 3-vector
     *  \param mu       Free direction (time and space) index (0-based)
     *  \param r_f      Final polarization index (0-based)
     *  \param r_i      Initial polarization index (0-based)
     */
    virtual MatElemRes_t<T1,T2> operator()(const std::string& FF_name,
					   const Real& mass_f, const ArrayInt& p_f, 
					   const Real& mass_i, const ArrayInt& p_i, 
					   const Array<int>& lorentz, 
					   int r_f, int r_i) const = 0;
  };


  //----------------------------------------------------------------------------------
  //! Params for matrix elements
  struct MatrixElementParams
  {
    MatrixElementParams();
    MatrixElementParams(XMLReader& in, const std::string& path);
    void writeXML(XMLWriter& in, const std::string& path) const;

    LatticeParam       lattice;            /*!< Holds lattice size and aniso*/
  };


  //! Reader
  void read(XMLReader& xml, const std::string& path, MatrixElementParams& param);

  //! Writer
  void write(XMLWriter& xml, const std::string& path, const MatrixElementParams& param);

} // namespace FF

#endif
