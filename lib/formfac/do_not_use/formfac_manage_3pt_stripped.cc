// $Id: formfac_manage_3pt_stripped.cc,v 2.0 2008/12/05 04:43:36 edwards Exp $
//
// Manage building-blocks


#include "formfac/formfac_manage_3pt_stripped.h"
#include "formfac/formfac_manage_factory.h"

namespace FF
{

  //----------------------------------------------------------------------------------
  //! Constructor
  Manage3PtFuncStripped::Manage3PtFuncStripped(ADAT::Handle<State3PtFunc> state_, 
					       const std::string& cache_file_,
					       int max_map_mb_) :
    Manage3PtFuncCache(cache_file_,max_map_mb_), state(state_) 
  {}


  //! Destructor
  Manage3PtFuncStripped::~Manage3PtFuncStripped() {}


  //! Read files
  void Manage3PtFuncStripped::do_read3pt(const ThreePtArg& arg)
  {
    std::cout << __func__ << ": entering" << std::endl;

    std::string filename = (*state)(0, arg);   // the cfg is ignored here
    std::cout << __func__ << ": read " << filename << std::endl;

    EnsemVectorComplexF thr;
    read(filename, thr);
    
    this->insert(arg, thr);

    std::cout << __func__ << ": exiting" << std::endl;
  }



#if 0
  //! Manager factory
  namespace FormfacManageEnv
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
	success &= TheManage3PtFuncFactory::Instance().registerObject(std::string("THREEPT_BB"),
								      matPi0VecCurPi0);

	registered = true;
      }

      return success;
    }

  }  // end namespace FormfacManageEnv
#endif

} // namespace FF
