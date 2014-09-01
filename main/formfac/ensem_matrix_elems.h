// -*- C++ -*-
// $Id: ensem_matrix_elems.h,v 2.0 2008/12/05 04:43:47 edwards Exp $
/*! \file
 * \brief Ensemble version of matrix elements
 */

#ifndef __ensem_matrix_elems_h__
#define __ensem_matrix_elems_h__

#include "adat/handle.h"
#include "matrix_elems.h"
#include "meson_matrix_elems.h"
#include "baryon_matrix_elems.h"

namespace FF
{
  using namespace ADAT;

  //----------------------------------------------------------------------------------
  //! Base class for Minkowski-space matrix elements
  template<typename T1, typename T2>
  class EnsemMatrixElement
  {
  public:  
    //! Virtual destructor
    virtual ~EnsemMatrixElement() {}

    //! The source state of this matrix element
    virtual std::string sourceName() const {return getMat().sourceName();}

    //! The sink state of this matrix element
    virtual std::string sinkName() const {return getMat().sinkName();}

    //! The number of source polarizations
    virtual int numSourcePolar() const {return getMat().numSourcePolar();}

    //! The number of sink polarizations
    virtual int numSinkPolar() const {return getMat().numSinkPolar();}

    //! The number of directions in the insertion
    virtual int numDir() const {return getMat().numDir();}

    //! The number of form-factors
    virtual int numFF() const {return getMat().numFF();}

    //! Return names for the form-factors
    virtual std::vector<std::string> getFFNames() const {return getMat().getFFNames();}

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
					   const EnsemReal& mass_f, const Array<int>& p_f, 
					   const EnsemReal& mass_i, const Array<int>& p_i, 
					   const Array<int>& lorentz, 
					   int r_f, int r_i) const = 0;
    
  protected:
    virtual const MatrixElement<typename EnsemScalar<T1>::Type_t,typename EnsemScalar<T2>::Type_t>& getMat() const = 0;
  };


  //----------------------------------------------------------------------------------
  //! Base class for meson-meson Minkowski-space matrix elements
  class EnsemMesonMatrixElement : public EnsemMatrixElement<EnsemComplex,EnsemComplex>
  {
  public:  
    // Typedefs to save typing
    typedef EnsemComplex  T1;
    typedef EnsemComplex  T2;

    //! Constructor
    EnsemMesonMatrixElement(Handle<MesonMatrixElement> mat_) : mat_elem(mat_) {}

    //! Destructor
    ~EnsemMesonMatrixElement() {}
    
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
					   const EnsemReal& mass_f, const Array<int>& p_f, 
					   const EnsemReal& mass_i, const Array<int>& p_i, 
					   const Array<int>& lorentz, 
					   int r_f, int r_i) const
      {
	MatElemRes_t<T1,T2> out;
	out.result.checkResize(__func__, mass_f.size(), mass_f.getEnsemType(), mass_i.size(), mass_i.getEnsemType());

	EnsemReal m_f(rescaleEnsemDown(mass_f));
	EnsemReal m_i(rescaleEnsemDown(mass_i));

	for(int bin=0; bin < out.result.size(); ++bin)
	{
	  Real ff = peekEnsem(m_f, bin);
	  Real ii = peekEnsem(m_i, bin);

	  MatElemRes_t<Complex,Complex> res = getMat()(FF_name, ff, p_f, ii, p_i, lorentz, r_f, r_i);

	  pokeEnsem(out.result, res.result, bin);
	}

	out.result = rescaleEnsemUp(out.result);
	return out;
      }
    
  protected:
    const MatrixElement<Complex,Complex>& getMat() const {return *mat_elem;}

  private:
    Handle<MesonMatrixElement> mat_elem;
  };


  //----------------------------------------------------------------------------------
  //! Base class for baryon-baryon Minkowski-space matrix elements
  class EnsemBaryonMatrixElement : public EnsemMatrixElement<EnsemSpinVector,EnsemSpinVector>
  {
  public:  
    // Typedefs to save typing
    typedef EnsemSpinVector  T1;
    typedef EnsemSpinVector  T2;

    //! Constructor
    EnsemBaryonMatrixElement(Handle<BaryonMatrixElement> mat_) : mat_elem(mat_) {}

    //! Destructor
    ~EnsemBaryonMatrixElement() {}
    
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
					   const EnsemReal& mass_f, const Array<int>& p_f, 
					   const EnsemReal& mass_i, const Array<int>& p_i, 
					   const Array<int>& lorentz, 
					   int r_f, int r_i) const
      {
	MatElemRes_t<T1,T2> out;
	out.result.checkResize(__func__, mass_f.size(), mass_f.getEnsemType(), mass_i.size(), mass_i.getEnsemType());
	out.source.checkResize(__func__, mass_f.size(), mass_f.getEnsemType(), mass_i.size(), mass_i.getEnsemType());
	out.sink.checkResize(__func__, mass_f.size(), mass_f.getEnsemType(), mass_i.size(), mass_i.getEnsemType());

	EnsemReal m_f(rescaleEnsemDown(mass_f));
	EnsemReal m_i(rescaleEnsemDown(mass_i));

	for(int bin=0; bin < out.result.size(); ++bin)
	{
	  Real ff = peekEnsem(m_f, bin);
	  Real ii = peekEnsem(m_i, bin);

	  MatElemRes_t<SpinVector,SpinVector> res = getMat()(FF_name, ff, p_f, ii, p_i, lorentz, r_f, r_i);
	  
	  pokeEnsem(out.result, res.result, bin);
	  pokeEnsem(out.source, res.source, bin);
	  pokeEnsem(out.sink, res.sink, bin);
	}

	out.result = rescaleEnsemUp(out.result);
	out.source = rescaleEnsemUp(out.source);
	out.sink   = rescaleEnsemUp(out.sink);
	return out;
      }
    
  protected:
    const MatrixElement<SpinVector,SpinVector>& getMat() const {return *mat_elem;}

  private:
    Handle<BaryonMatrixElement> mat_elem;
  };



} // namespace FF

#endif
