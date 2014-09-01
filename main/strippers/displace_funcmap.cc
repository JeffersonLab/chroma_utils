// $Id: displace_funcmap.cc,v 2.0 2008/12/05 04:44:03 edwards Exp $
/*! \file
 *  \brief Displacement function map
 */

#include "strippers.h"
#include "displace_funcmap.h"
#include <string>
#include "io/adat_xmlio.h"

namespace Strippers
{
  using namespace Util;

  //! Displacement functions
  namespace DisplacementCallMapEnv
  { 
    //! Anonymous namespace
    namespace
    {
      const std::string disp_path = "/Displacement";

      //! Params for derivative quark displacement
      struct Params
      {
	Params() {deriv_length = 0;}
	Params(XMLReader& xml, const std::string& path)
	  {
	    XMLReader paramtop(xml, path);

	    int version;
	    read(paramtop, "version", version);

	    switch (version) 
	    {
	    case 1:
	      break;

	    default:
	      std::cerr << __func__ << ": parameter version " << version 
			<< " unsupported." << std::endl;
	      exit(1);
	    }

	    read(paramtop, "DisplacementType", displacement_type);
	    read(paramtop, "deriv_length", deriv_length);
	  }

	std::string      displacement_type;    /*!< displacement type */

	int              deriv_length;         /*!< Displacement length in derivative */
      };

      //! Deriv meson source parameters
      struct ParamsDir
      {
	ParamsDir() {deriv_dir = deriv_length = 0;}
	ParamsDir(XMLReader& xml, const std::string& path)
	  {
	    XMLReader paramtop(xml, path);

	    int version;
	    read(paramtop, "version", version);

	    switch (version) 
	    {
	    case 1:
	      break;

	    default:
	      std::cerr << __func__ << ": parameter version " << version 
			<< " unsupported." << std::endl;
	      exit(1);
	    }

	    read(paramtop, "DisplacementType", displacement_type);
	    read(paramtop, "deriv_dir", deriv_dir);
	    read(paramtop, "deriv_length", deriv_length);
	  }

	std::string      displacement_type;    /*!< displacement type */

	int              deriv_dir;            /*!< Polarization direction */
	int              deriv_length;         /*!< Displacement length in derivative */
      };


      // Read parameters
      void read(XMLReader& xml, const std::string& path, Params& param)
      {
	Params tmp(xml, path);
	param = tmp;
      }

      // Read parameters
      void read(XMLReader& xml, const std::string& path, ParamsDir& param)
      {
	ParamsDir tmp(xml, path);
	param = tmp;
      }


      //! Map an int to a direction for a T1 or T2 rep
      std::string TtoDir[3] = {"_x", "_y", "_z"};

      //! Map an int to a direction for an E rep
      std::string EtoDir[2] = {"_0", "_1"};


      DisplacementNames_t mesPionxNablaT1Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "1pm" + TtoDir[params.deriv_dir];
	names.wavetype = "pixNablaT1";
	return names;
      }

      DisplacementNames_t mesA0xNablaT1Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "1mm" + TtoDir[params.deriv_dir];
	names.wavetype = "a0xNablaT1";
	return names;
      }

      DisplacementNames_t mesA02xNablaT1Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "1mp" + TtoDir[params.deriv_dir];
	names.wavetype = "a02xNablaT1";
	return names;
      }

      DisplacementNames_t mesRhoxNablaA1Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "0pp";
	names.wavetype = "rhoxNablaA1";
	return names;
      }

      DisplacementNames_t mesRhoxNablaT1Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "1pp" + TtoDir[params.deriv_dir];
	names.wavetype = "rhoxNablaT1";
	return names;
      }

      DisplacementNames_t mesRhoxNablaT2Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "2pp" + TtoDir[params.deriv_dir];
	names.wavetype = "rhoxNablaT2";
	return names;
      }

      DisplacementNames_t mesA1xNablaA1Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "0mm";
	names.wavetype = "a1xNablaA1";
	return names;
      }

      DisplacementNames_t mesA1xNablaT2Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "2mm" + TtoDir[params.deriv_dir];
	names.wavetype = "a1xNablaT2";
	return names;
      }

      DisplacementNames_t mesA1xNablaEDisplace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "2mm" + EtoDir[params.deriv_dir];
	names.wavetype = "a1xNablaE";
	return names;
      }

      DisplacementNames_t mesB1xNablaT1Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "1mp" + TtoDir[params.deriv_dir];
	names.wavetype = "b1xNablaT1";
	return names;
      }

      DisplacementNames_t mesA02xDT2Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "2pm" + TtoDir[params.deriv_dir];
	names.wavetype = "a02xDT2";
	return names;
      }

      DisplacementNames_t mesA1xDA2Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "3pp";
	names.wavetype = "a1xDA2";
	return names;
      }

      DisplacementNames_t mesA1xDEDisplace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "2pp" + EtoDir[params.deriv_dir];
	names.wavetype = "a1xDE";
	return names;
      }

      DisplacementNames_t mesA1xDT1Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "1pp" + TtoDir[params.deriv_dir];
	names.wavetype = "a1xDT1";
	return names;
      }

      DisplacementNames_t mesA1xDT2Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "2pp" + TtoDir[params.deriv_dir];
	names.wavetype = "a1xDT2";
	return names;
      }

      DisplacementNames_t mesB1xDA2Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "3pm";
	names.wavetype = "b1xDA2";
	return names;
      }

      DisplacementNames_t mesB1xDEDisplace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "2pm" + EtoDir[params.deriv_dir];
	names.wavetype = "b1xDE";
	return names;
      }

      DisplacementNames_t mesB1xDT1Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "1pm" + TtoDir[params.deriv_dir];
	names.wavetype = "b1xDT1";
	return names;
      }

      DisplacementNames_t mesB1xDT2Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "3pm" + TtoDir[params.deriv_dir];
	names.wavetype = "b1xDT2";
	return names;
      }

      DisplacementNames_t mesRhoxDA2Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "3mm";
	names.wavetype = "rhoxDA2";
	return names;
      }

      DisplacementNames_t mesRhoxDT1Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "1mm" + TtoDir[params.deriv_dir];
	names.wavetype = "rhoxDT1";
	return names;
      }

      DisplacementNames_t mesRhoxDT2Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "2mm" + TtoDir[params.deriv_dir];
	names.wavetype = "rhoxDT2";
	return names;
      }

      DisplacementNames_t mesPionxDT2Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "2mp" + TtoDir[params.deriv_dir];
	names.wavetype = "pixDT2";
	return names;
      }

      DisplacementNames_t mesPionxBT1Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "1mm" + TtoDir[params.deriv_dir];
	names.wavetype = "pixBT1";
	return names;
      }

      DisplacementNames_t mesRhoxBT1Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "1mp" + TtoDir[params.deriv_dir];
	names.wavetype = "rhoxBT1";
	return names;
      }

      DisplacementNames_t mesRhoxBT2Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "2mp" + TtoDir[params.deriv_dir];
	names.wavetype = "rhoxBT2";
	return names;
      }

      DisplacementNames_t mesA1xBA1Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "0pm";
	names.wavetype = "a1xBA1";
	return names;
      }

      DisplacementNames_t mesA1xBT1Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "1pm" + TtoDir[params.deriv_dir];
	names.wavetype = "a1xBT1";
	return names;
      }

      DisplacementNames_t mesA1xBT2Displace(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	ParamsDir params(displacetop, disp_path);

	names.particle = "2pm" + TtoDir[params.deriv_dir];
	names.wavetype = "a1xBT2";
	return names;
      }

    } // end anonymous namespace


    //! Another namespace for gamma insertions
    namespace 
    {
      //! Params for gamma insertions
      struct GammaParams
      {
	GammaParams() {gamma = 0;}
	GammaParams(XMLReader& xml, const std::string& path)
	  {
	    XMLReader paramtop(xml, path);

	    int version;
	    read(paramtop, "version", version);

	    switch (version) 
	    {
	    case 1:
	      break;

	    default:
	      std::cerr << __func__ << ": parameter version " << version 
			<< " unsupported." << std::endl;
	      exit(1);
	    }

	    read(paramtop, "DisplacementType", displacement_type);
	    read(paramtop, "gamma", gamma);
	  }

	std::string      displacement_type;    /*!< displacement type */

	int              gamma;                /*!< Gamma insertion */
      };


      //! Map an int to a direction for an E rep
      std::string GammaToName[16] = {"a0", "1mm_x", "1mm_y", "b1_z", "1mm_z", "b1_y", "b1_x", "pion",
      "a0", "1mm_x", "1mm_y", "a1_z", "1mm_z", "a1_y", "a1_x", "pion"};

      //! Map an int to a direction for an E rep
      std::string MesonState[16] = {"1", "1", "1", "1", "1", "1", "1", "2", 
      "2", "2", "2", "1", "2", "1", "1", "1"};


      DisplacementNames_t mesGammaInsertion(const std::string& displacement)
      {
	DisplacementNames_t names;

	std::stringstream  is_d(displacement);
	XMLReader  displacetop(is_d);
	GammaParams params(displacetop, disp_path);

	names.particle = GammaToName[params.gamma];
	names.wavetype = MesonState[params.gamma];
	return names;
      }


    }


    bool registerAll(void) 
    {
      bool success = true;
      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("PIONxNABLA_T1-DERIV"),
								     mesPionxNablaT1Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("A0xNABLA_T1-DERIV"),
								     mesA0xNablaT1Displace);
      
      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("A0_2xNABLA_T1-DERIV"),
								     mesA02xNablaT1Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("RHOxNABLA_A1-DERIV"),
								     mesRhoxNablaA1Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("RHOxNABLA_T1-DERIV"),
								     mesRhoxNablaT1Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("RHOxNABLA_T2-DERIV"),
								     mesRhoxNablaT2Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("A1xNABLA_A1-DERIV"),
								     mesA1xNablaA1Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("A1xNABLA_T2-DERIV"),
								     mesA1xNablaT2Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("A1xNABLA_E-DERIV"),
								     mesA1xNablaEDisplace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("B1xNABLA_T1-DERIV"),
								     mesB1xNablaT1Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("A0_2xD_T2-DERIV"),
								     mesA02xDT2Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("A1xD_A2-DERIV"),
								     mesA1xDA2Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("A1xD_E-DERIV"),
								     mesA1xDEDisplace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("A1xD_T1-DERIV"),
								     mesA1xDT1Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("A1xD_T2-DERIV"),
								     mesA1xDT2Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("B1xD_A2-DERIV"),
								     mesB1xDA2Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("B1xD_E-DERIV"),
								     mesB1xDEDisplace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("B1xD_T1-DERIV"),
								     mesB1xDT1Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("B1xD_T2-DERIV"),
								     mesB1xDT2Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("RHOxD_A2-DERIV"),
								     mesRhoxDA2Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("RHOxD_T1-DERIV"),
								     mesRhoxDT1Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("RHOxD_T2-DERIV"),
								     mesRhoxDT2Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("PIONxD_T2-DERIV"),
								     mesPionxDT2Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("PIONxB_T1-DERIV"),
								     mesPionxBT1Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("RHOxB_T1-DERIV"),
								     mesRhoxBT1Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("RHOxB_T2-DERIV"),
								     mesRhoxBT2Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("A1xB_A1-DERIV"),
								     mesA1xBA1Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("A1xB_T1-DERIV"),
								     mesA1xBT1Displace);

      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("A1xB_T2-DERIV"),
								     mesA1xBT2Displace);

      // Gamma insertions
      success &= TheDisplacementFuncMap::Instance().registerFunction(std::string("GAMMA_INSERTION"),
								     mesGammaInsertion);

      return success;
    }

    bool registered = registerAll();
  } // namespace DisplacementCallMapEnv

}

