// $Id: meson_matrix_elems.cc,v 2.0 2008/12/05 04:43:48 edwards Exp $
/*! \file
 * \brief Matrix elements
 */

#include "wavefuncs.h"
#include "matrix_elems.h"
#include "matrix_elems_factory.h"
#include <iostream>

#include <math.h>

using namespace std;

namespace FF
{
  //---------------------------------------------------------------------------
  // Return names for the form-factors
  std::vector<std::string> 
  Pi0VecCurPi0MatElem::getFFNames() const
  {
    std::vector<string> names;
    names.push_back("F_pi");

    return names;
  }

  // The matrix element value
  MatElemRes_t<Complex,Complex>
  Pi0VecCurPi0MatElem::operator()(const std::string& FF_name,
				  const Real& mass_f, const ArrayInt& p_f, 
				  const Real& mass_i, const ArrayInt& p_i, 
				  const Array<int>& lorentz, 
				  int r_f, int r_i) const
  {
    MatElemRes_t<Complex,Complex> value;

    if (lorentz.size() != 1)
    {
      std::cerr << getFFNames()[0] << ": unexpected lorentz index size= " << lorentz.size() << std::endl;
      exit(1);
    }
    int mu = lorentz[0];

    Array<Real> pp_f = make4Vec(mass_f,p_f,params.lattice);
    Array<Real> pp_i = make4Vec(mass_i,p_i,params.lattice);
    Real tmp = pp_f[mu] + pp_i[mu];

    if (FF_name == "F_pi")
    {
      value.result = tmp;
    }
    else
    {
      value.result = zero;
    }

    return value; 
  }


 
  //! Meson sources
  namespace MesonMatrixElementEnv
  { 
    //! Anonymous namespace
    namespace
    {
      bool registered = false;

      //-------------------- callback functions ---------------------------------------

      //! pi0|V|pi0
      MesonMatrixElement* matPi0VecCurPi0(XMLReader& xml_in,
					  const std::string& path)
      {
	return new Pi0VecCurPi0MatElem(MatrixElementParams(xml_in, path));
      }

    }  // anonymous namespace


    // Register all the possible matrix elements
    bool registerAll(void) 
    {
      bool success = true;

      if (! registered)
      {
	//! Register all the factories
	success &= TheMesonMatrixElementFactory::Instance().registerObject(std::string("pi0|VecCur|pi0"),
									   matPi0VecCurPi0);

	registered = true;
      }

      return success;
    }

  }  // end namespace MatrixElementEnv

} // namespace FF
