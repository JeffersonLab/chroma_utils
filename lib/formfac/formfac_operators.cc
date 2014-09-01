// $Id: formfac_operators.cc,v 2.0 2008/12/05 04:43:36 edwards Exp $
/*! \file
 * \brief Operators and wavefunction overlaps
 */

#include "operators.h"
#include "operators_factory.h"
#include "formfac/formfac.h"

namespace FF
{
  using namespace Util;

  // Read parameters
  void read(XMLReader& xml, const std::string& path, OperatorParams& param)
  {
    OperatorParams tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const std::string& path, const OperatorParams& param)
  {
    param.writeXML(xml, path);
  }

  //! Initialize
  OperatorParams::OperatorParams()
  {
  }


  //! Read parameters
  OperatorParams::OperatorParams(XMLReader& xml, const std::string& path)
  {
    XMLReader paramtop(xml, path);

//      read(paramtop, "OperatorType",  operator_type);
    read(paramtop, "LatticeParam", lattice);
  }


  // Writer
  void OperatorParams::writeXML(XMLWriter& xml, const std::string& path) const
  {
    push(xml, path);

//      write(xml, "OperatorType", operator_type);
    write(xml, "LatticeParam", lattice);

    pop(xml);
  }



  //! Anonymous namespace
  namespace
  {
    //-------------------- callback functions ---------------------------------------

    //----------------------------
    // Gamma(n) variants

    //! Construct (Gamma(0)) operator
    Operator<Complex>* mesGamma0Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesGamma0Operator(OperatorParams(xml_in, path));
    }

    //! Construct (Gamma(1)) operator
    Operator<Complex>* mesGamma1Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesGamma1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (Gamma(2)) operator
    Operator<Complex>* mesGamma2Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesGamma2Operator(OperatorParams(xml_in, path));
    }

    //! Construct (Gamma(3)) operator
    Operator<Complex>* mesGamma3Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesGamma3Operator(OperatorParams(xml_in, path));
    }

    //! Construct (Gamma(4)) operator
    Operator<Complex>* mesGamma4Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesGamma4Operator(OperatorParams(xml_in, path));
    }

    //! Construct (Gamma(5)) operator
    Operator<Complex>* mesGamma5Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesGamma5Operator(OperatorParams(xml_in, path));
    }

    //! Construct (Gamma(6)) operator
    Operator<Complex>* mesGamma6Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesGamma6Operator(OperatorParams(xml_in, path));
    }

    //! Construct (Gamma(7)) operator
    Operator<Complex>* mesGamma7Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesGamma7Operator(OperatorParams(xml_in, path));
    }

    //! Construct (Gamma(8)) operator
    Operator<Complex>* mesGamma8Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesGamma8Operator(OperatorParams(xml_in, path));
    }

    //! Construct (Gamma(9)) operator
    Operator<Complex>* mesGamma9Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesGamma9Operator(OperatorParams(xml_in, path));
    }

    //! Construct (Gamma(10)) operator
    Operator<Complex>* mesGamma10Operator(XMLReader& xml_in,
					  const std::string& path)
    {
      return new MesGamma10Operator(OperatorParams(xml_in, path));
    }

    //! Construct (Gamma(11)) operator
    Operator<Complex>* mesGamma11Operator(XMLReader& xml_in,
					  const std::string& path)
    {
      return new MesGamma11Operator(OperatorParams(xml_in, path));
    }

    //! Construct (Gamma(12)) operator
    Operator<Complex>* mesGamma12Operator(XMLReader& xml_in,
					  const std::string& path)
    {
      return new MesGamma12Operator(OperatorParams(xml_in, path));
    }

    //! Construct (Gamma(13)) operator
    Operator<Complex>* mesGamma13Operator(XMLReader& xml_in,
					  const std::string& path)
    {
      return new MesGamma13Operator(OperatorParams(xml_in, path));
    }

    //! Construct (Gamma(14)) operator
    Operator<Complex>* mesGamma14Operator(XMLReader& xml_in,
					  const std::string& path)
    {
      return new MesGamma14Operator(OperatorParams(xml_in, path));
    }

    //! Construct (Gamma(15)) operator
    Operator<Complex>* mesGamma15Operator(XMLReader& xml_in,
					  const std::string& path)
    {
      return new MesGamma15Operator(OperatorParams(xml_in, path));
    }


    //----------------------------
    // Simple gamma's by name

    //! Construct (a0) operator
    Operator<Complex>* mesA0A1Operator(XMLReader& xml_in,
				       const std::string& path)
    {
      return new MesA0A1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (pi_0) operator
    Operator<Complex>* mesPionA1Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesPionA1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (rho_1) operator
    Operator<Complex>* mesRhoT1Operator(XMLReader& xml_in,
					const std::string& path)
    {
      return new MesRhoT1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (a_1) operator
    Operator<Complex>* mesA1T1Operator(XMLReader& xml_in,
				       const std::string& path)
    {
      return new MesA1T1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (b_1) operator
    Operator<Complex>* mesB1T1Operator(XMLReader& xml_in,
				       const std::string& path)
    {
      return new MesB1T1Operator(OperatorParams(xml_in, path));
    }


    //-------------------------------
    //! Construct (PionxNabla_T1) operator
    Operator<Complex>* mesPionxNablaT1Operator(XMLReader& xml_in,
					       const std::string& path)
    {
      return new MesPionxNablaT1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (A0xNabla_T1) operator
    Operator<Complex>* mesA0xNablaT1Operator(XMLReader& xml_in,
					     const std::string& path)
    {
      return new MesA0xNablaT1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (A0_2xNabla_T1) operator
    Operator<Complex>* mesA02xNablaT1Operator(XMLReader& xml_in,
					      const std::string& path)
    {
      return new MesA02xNablaT1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (RhoxNabla_A1) operator
    Operator<Complex>* mesRhoxNablaA1Operator(XMLReader& xml_in,
					      const std::string& path)
    {
      return new MesRhoxNablaA1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (RhoxNabla_T1) operator
    Operator<Complex>* mesRhoxNablaT1Operator(XMLReader& xml_in,
					      const std::string& path)
    {
      return new MesRhoxNablaT1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (RhoxNabla_T2) operator
    Operator<Complex>* mesRhoxNablaT2Operator(XMLReader& xml_in,
					      const std::string& path)
    {
      return new MesRhoxNablaT2Operator(OperatorParams(xml_in, path));
    }

    //! Construct (RhoxNabla_E) operator
    Operator<Complex>* mesRhoxNablaEOperator(XMLReader& xml_in,
					     const std::string& path)
    {
      return new MesRhoxNablaEOperator(OperatorParams(xml_in, path));
    }

    //! Construct (A1xNabla_A1) operator
    Operator<Complex>* mesA1xNablaA1Operator(XMLReader& xml_in,
					     const std::string& path)
    {
      return new MesA1xNablaA1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (A1xNabla_T1) operator
    Operator<Complex>* mesA1xNablaT1Operator(XMLReader& xml_in,
					     const std::string& path)
    {
      return new MesA1xNablaT1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (A1xNabla_T2) operator
    Operator<Complex>* mesA1xNablaT2Operator(XMLReader& xml_in,
					     const std::string& path)
    {
      return new MesA1xNablaT2Operator(OperatorParams(xml_in, path));
    }

    //! Construct (A1xNabla_E) operator
    Operator<Complex>* mesA1xNablaEOperator(XMLReader& xml_in,
					    const std::string& path)
    {
      return new MesA1xNablaEOperator(OperatorParams(xml_in, path));
    }

    //! Construct (B1xNabla_A1) operator
    Operator<Complex>* mesB1xNablaA1Operator(XMLReader& xml_in,
					     const std::string& path)
    {
      return new MesB1xNablaA1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (B1xNabla_T1) operator
    Operator<Complex>* mesB1xNablaT1Operator(XMLReader& xml_in,
					     const std::string& path)
    {
      return new MesB1xNablaT1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (B1xNabla_T2) operator
    Operator<Complex>* mesB1xNablaT2Operator(XMLReader& xml_in,
					     const std::string& path)
    {
      return new MesB1xNablaT2Operator(OperatorParams(xml_in, path));
    }

    //! Construct (B1xNabla_E) operator
    Operator<Complex>* mesB1xNablaEOperator(XMLReader& xml_in,
					    const std::string& path)
    {
      return new MesB1xNablaEOperator(OperatorParams(xml_in, path));
    }

    //! Construct (A0_2xD_T2) operator
    Operator<Complex>* mesA02xDT2Operator(XMLReader& xml_in,
					  const std::string& path)
    {
      return new MesA02xDT2Operator(OperatorParams(xml_in, path));
    }

    //! Construct (A1xD_A2) operator
    Operator<Complex>* mesA1xDA2Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesA1xDA2Operator(OperatorParams(xml_in, path));
    }

    //! Construct (A1xD_E) operator
    Operator<Complex>* mesA1xDEOperator(XMLReader& xml_in,
					const std::string& path)
    {
      return new MesA1xDEOperator(OperatorParams(xml_in, path));
    }

    //! Construct (A1xD_T1) operator
    Operator<Complex>* mesA1xDT1Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesA1xDT1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (A1xD_T2) operator
    Operator<Complex>* mesA1xDT2Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesA1xDT2Operator(OperatorParams(xml_in, path));
    }


    //! Construct (B1xD_A2) operator
    Operator<Complex>* mesB1xDA2Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesB1xDA2Operator(OperatorParams(xml_in, path));
    }

    //! Construct (B1xD_E) operator
    Operator<Complex>* mesB1xDEOperator(XMLReader& xml_in,
					const std::string& path)
    {
      return new MesB1xDEOperator(OperatorParams(xml_in, path));
    }

    //! Construct (B1xD_T1) operator
    Operator<Complex>* mesB1xDT1Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesB1xDT1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (B1xD_T2) operator
    Operator<Complex>* mesB1xDT2Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesB1xDT2Operator(OperatorParams(xml_in, path));
    }

    //! Construct (RhoxD_A2) operator
    Operator<Complex>* mesRhoxDA2Operator(XMLReader& xml_in,
					  const std::string& path)
    {
      return new MesRhoxDA2Operator(OperatorParams(xml_in, path));
    }

    //! Construct (RhoxD_T1) operator
    Operator<Complex>* mesRhoxDT1Operator(XMLReader& xml_in,
					  const std::string& path)
    {
      return new MesRhoxDT1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (RhoxD_T2) operator
    Operator<Complex>* mesRhoxDT2Operator(XMLReader& xml_in,
					  const std::string& path)
    {
      return new MesRhoxDT2Operator(OperatorParams(xml_in, path));
    }

    //! Construct (RhoxD_E) operator
    Operator<Complex>* mesRhoxDEOperator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesRhoxDEOperator(OperatorParams(xml_in, path));
    }

    //! Construct (PionxD_T2) operator
    Operator<Complex>* mesPionxDT2Operator(XMLReader& xml_in,
					   const std::string& path)
    {
      return new MesPionxDT2Operator(OperatorParams(xml_in, path));
    }

    //! Construct (Pion2xD_T2) operator
    Operator<Complex>* mesPion2xDT2Operator(XMLReader& xml_in,
					    const std::string& path)
    {
      return new MesPion2xDT2Operator(OperatorParams(xml_in, path));
    }

    //! Construct (A0xD_T2) operator
    Operator<Complex>* mesA0xDT2Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesA0xDT2Operator(OperatorParams(xml_in, path));
    }

    //! Construct (A0xB_T1) operator
    Operator<Complex>* mesA0xBT1Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesA0xBT1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (A02xB_T1) operator
    Operator<Complex>* mesA02xBT1Operator(XMLReader& xml_in,
					  const std::string& path)
    {
      return new MesA02xBT1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (PionxB_T1) operator
    Operator<Complex>* mesPionxBT1Operator(XMLReader& xml_in,
					   const std::string& path)
    {
      return new MesPionxBT1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (RhoxB_A1) operator
    Operator<Complex>* mesRhoxBA1Operator(XMLReader& xml_in,
					  const std::string& path)
    {
      return new MesRhoxBA1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (RhoxB_T1) operator
    Operator<Complex>* mesRhoxBT1Operator(XMLReader& xml_in,
					  const std::string& path)
    {
      return new MesRhoxBT1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (RhoxB_T2) operator
    Operator<Complex>* mesRhoxBT2Operator(XMLReader& xml_in,
					  const std::string& path)
    {
      return new MesRhoxBT2Operator(OperatorParams(xml_in, path));
    }

    //! Construct (RhoxB_E) operator
    Operator<Complex>* mesRhoxBEOperator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesRhoxBEOperator(OperatorParams(xml_in, path));
    }


    //! Construct (A1xB_A1) operator
    Operator<Complex>* mesA1xBA1Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesA1xBA1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (A1xB_T1) operator
    Operator<Complex>* mesA1xBT1Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesA1xBT1Operator(OperatorParams(xml_in, path));
    }

    //! Construct (A1xB_T2) operator
    Operator<Complex>* mesA1xBT2Operator(XMLReader& xml_in,
					 const std::string& path)
    {
      return new MesA1xBT2Operator(OperatorParams(xml_in, path));
    }

    //! Construct (A1xB_E) operator
    Operator<Complex>* mesA1xBEOperator(XMLReader& xml_in,
					const std::string& path)
    {
      return new MesA1xBEOperator(OperatorParams(xml_in, path));
    }
  } // end anonymous namespace


    //----------------------------------------------------------------------------------
    // Simple \Gamma(n) operators

    // Construct (Gamma(0)) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesGamma0Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesGamma0A0WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (Gamma(1)) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesGamma1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesGamma1B0WaveFunc(wvfParams));
    wvf.push_back(new MesGamma1Rho1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (Gamma(2)) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesGamma2Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesGamma2B0WaveFunc(wvfParams));
    wvf.push_back(new MesGamma2Rho1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (Gamma(3)) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesGamma3Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesGamma3B1WaveFunc(wvfParams));
    wvf.push_back(new MesGamma3Rho1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (Gamma(4)) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesGamma4Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesGamma4B0WaveFunc(wvfParams));
    wvf.push_back(new MesGamma4Rho1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (Gamma(5)) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesGamma5Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesGamma5B1WaveFunc(wvfParams));
    wvf.push_back(new MesGamma5Rho1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (Gamma(6)) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesGamma6Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesGamma6B1WaveFunc(wvfParams));
    wvf.push_back(new MesGamma6Rho1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (Gamma(7)) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesGamma7Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesGamma7Pi0WaveFunc(wvfParams));
    wvf.push_back(new MesGamma7A1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (Gamma(8)) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesGamma8Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesGamma8B0WaveFunc(wvfParams));
    wvf.push_back(new MesGamma8Rho1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (Gamma(9)) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesGamma9Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesGamma9B1WaveFunc(wvfParams));
    wvf.push_back(new MesGamma9Rho1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (Gamma(10)) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesGamma10Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesGamma10B1WaveFunc(wvfParams));
    wvf.push_back(new MesGamma10Rho1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (Gamma(11)) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesGamma11Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesGamma11Pi0WaveFunc(wvfParams));
    wvf.push_back(new MesGamma11A1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (Gamma(12)) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesGamma12Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesGamma12B1WaveFunc(wvfParams));
    wvf.push_back(new MesGamma12Rho1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (Gamma(13)) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesGamma13Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesGamma13Pi0WaveFunc(wvfParams));
    wvf.push_back(new MesGamma13A1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (Gamma(14)) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesGamma14Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesGamma14Pi0WaveFunc(wvfParams));
    wvf.push_back(new MesGamma14A1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (Gamma(15)) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesGamma15Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesGamma15Pi0WaveFunc(wvfParams));

    return wvf;
  }


  //----------------------------------------------------------------------------------
  // Construct (a_0) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA0A1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA0A1A0WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (pi_0) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesPionA1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesPionA1Pi0WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (rho_1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesRhoT1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesRhoT1B0WaveFunc(wvfParams));
    wvf.push_back(new MesRhoT1Rho1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (a_1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA1T1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA1T1Pi0WaveFunc(wvfParams));
    wvf.push_back(new MesA1T1A1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (b_1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesB1T1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesB1T1B1WaveFunc(wvfParams));
    wvf.push_back(new MesB1T1Rho1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (PionxNabla_T1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesPionxNablaT1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesPionxNablaT1Rho0WaveFunc(wvfParams));
    wvf.push_back(new MesPionxNablaT1B1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (A0xNabla_T1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA0xNablaT1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA0xNablaT1B0WaveFunc(wvfParams));
    wvf.push_back(new MesA0xNablaT1Rho1WaveFunc(wvfParams));

    return wvf;
      
  }


  // Construct (A0_2xNabla_T1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA02xNablaT1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA02xNablaT1A0WaveFunc(wvfParams));
    wvf.push_back(new MesA02xNablaT1A1WaveFunc(wvfParams));
    wvf.push_back(new MesA02xNablaT1A2WaveFunc(wvfParams));
    wvf.push_back(new MesA02xNablaT1Pi1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (RhoxNabla_A1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesRhoxNablaA1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesRhoxNablaA1A0WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxNablaA1A2WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxNablaA1Pi1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (RhoxNabla_T1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesRhoxNablaT1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesRhoxNablaT1A1WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxNablaT1Pi1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (RhoxNabla_T2) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesRhoxNablaT2Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesRhoxNablaT2A0WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxNablaT2A2WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxNablaT2Pi1WaveFunc(wvfParams));

    return wvf;
  }

  // Construct (RhoxNabla_E) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesRhoxNablaEOperator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesRhoxNablaEA0WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxNablaEA2WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxNablaEPi1WaveFunc(wvfParams));

    return wvf;
  }

  // Construct (A1xNabla_A1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA1xNablaA1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA1xNablaA1Rho0WaveFunc(wvfParams));
    wvf.push_back(new MesA1xNablaA1Rho2WaveFunc(wvfParams));
    wvf.push_back(new MesA1xNablaA1B1WaveFunc(wvfParams));

    return wvf;
  }

  // Construct (A1xNabla_T1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA1xNablaT1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA1xNablaT1Rho1WaveFunc(wvfParams));
    wvf.push_back(new MesA1xNablaT1B1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (A1xNabla_T2) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA1xNablaT2Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA1xNablaT2Rho0WaveFunc(wvfParams));
    wvf.push_back(new MesA1xNablaT2Rho2WaveFunc(wvfParams));
    wvf.push_back(new MesA1xNablaT2B1WaveFunc(wvfParams));

    return wvf;
  }


  // Construct (A1xNabla_E) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA1xNablaEOperator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA1xNablaERho0WaveFunc(wvfParams));
    wvf.push_back(new MesA1xNablaERho2WaveFunc(wvfParams));
    wvf.push_back(new MesA1xNablaEB1WaveFunc(wvfParams));

    return wvf;
  }

  // Construct (B1xNabla_A1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesB1xNablaA1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesB1xNablaA1Pi0WaveFunc(wvfParams));
    wvf.push_back(new MesB1xNablaA1A1WaveFunc(wvfParams));
    wvf.push_back(new MesB1xNablaA1Pi2WaveFunc(wvfParams));

    return wvf;
  }

  // Construct (B1xNabla_T1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesB1xNablaT1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesB1xNablaT1A0WaveFunc(wvfParams));
    wvf.push_back(new MesB1xNablaT1A1WaveFunc(wvfParams));
    wvf.push_back(new MesB1xNablaT1Pi1WaveFunc(wvfParams));
    wvf.push_back(new MesB1xNablaT1A2WaveFunc(wvfParams));
    wvf.push_back(new MesB1xNablaT1Pi2WaveFunc(wvfParams));

    return wvf;
  }

  // Construct (B1xNabla_T2) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesB1xNablaT2Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesB1xNablaT2A1WaveFunc(wvfParams));
    wvf.push_back(new MesB1xNablaT2Pi1WaveFunc(wvfParams));
    wvf.push_back(new MesB1xNablaT2A2WaveFunc(wvfParams));
    wvf.push_back(new MesB1xNablaT2Pi2WaveFunc(wvfParams));

    return wvf;
  }

  // Construct (B1xNabla_E) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesB1xNablaEOperator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesB1xNablaEA1WaveFunc(wvfParams));
    wvf.push_back(new MesB1xNablaEPi1WaveFunc(wvfParams));
    wvf.push_back(new MesB1xNablaEA2WaveFunc(wvfParams));
    wvf.push_back(new MesB1xNablaEPi2WaveFunc(wvfParams));

    return wvf;
  }

  // Construct (A0_2xD_T2) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA02xDT2Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA02xDT2B0WaveFunc(wvfParams));
    wvf.push_back(new MesA02xDT2B1WaveFunc(wvfParams));
    wvf.push_back(new MesA02xDT2Rho1WaveFunc(wvfParams));
    wvf.push_back(new MesA02xDT2B2WaveFunc(wvfParams));
    wvf.push_back(new MesA02xDT2Rho2WaveFunc(wvfParams));
    wvf.push_back(new MesA02xDT2Rho3WaveFunc(wvfParams));
    return wvf;
  }


  // Construct (A1xD_A2) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA1xDA2Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA1xDA2Pi2WaveFunc(wvfParams));
    wvf.push_back(new MesA1xDA2A3WaveFunc(wvfParams));

    return wvf;
  }

  // Construct (A1xD_T1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA1xDT1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA1xDT1Pi0WaveFunc(wvfParams));
    wvf.push_back(new MesA1xDT1Pi1WaveFunc(wvfParams));
    wvf.push_back(new MesA1xDT1A1WaveFunc(wvfParams));
    wvf.push_back(new MesA1xDT1Pi2WaveFunc(wvfParams));
    wvf.push_back(new MesA1xDT1A2WaveFunc(wvfParams));
    wvf.push_back(new MesA1xDT1A3WaveFunc(wvfParams));


    return wvf;
  }

  // Construct (A1xD_T2) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA1xDT2Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA1xDT2Pi0WaveFunc(wvfParams));
    wvf.push_back(new MesA1xDT2Pi1WaveFunc(wvfParams));
    wvf.push_back(new MesA1xDT2A1WaveFunc(wvfParams));
    wvf.push_back(new MesA1xDT2Pi2WaveFunc(wvfParams));
    wvf.push_back(new MesA1xDT2A2WaveFunc(wvfParams));
    wvf.push_back(new MesA1xDT2A3WaveFunc(wvfParams));


    return wvf;
  }

  // Construct (A1xD_E) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA1xDEOperator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA1xDEPi1WaveFunc(wvfParams));
    wvf.push_back(new MesA1xDEA1WaveFunc(wvfParams));
    wvf.push_back(new MesA1xDEPi2WaveFunc(wvfParams));
    wvf.push_back(new MesA1xDEA2WaveFunc(wvfParams));

    return wvf;
  }

  // Construct (B1xD_A2) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesB1xDA2Operator::overlaps() const
  {
  }


  // Construct (B1xD_E) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesB1xDEOperator::overlaps() const
  {
  }


  // Construct (B1xD_T1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesB1xDT1Operator::overlaps() const
  {
  }


  // Construct (B1xD_T2) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesB1xDT2Operator::overlaps() const
  {
  }


  // Construct (RhoxD_A2) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesRhoxDA2Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesRhoxDA2B2WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxDA2Rho3WaveFunc(wvfParams));
    return wvf;
  }


  //! Construct (RhoxD_T1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesRhoxDT1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesRhoxDT1B0WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxDT1B1WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxDT1Rho1WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxDT1B2WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxDT1Rho2WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxDT1Rho3WaveFunc(wvfParams));

    return wvf;
  }


  //! Construct (RhoxD_T2) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesRhoxDT2Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesRhoxDT2B0WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxDT2B1WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxDT2Rho1WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxDT2B2WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxDT2Rho2WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxDT2Rho3WaveFunc(wvfParams));

    return wvf;
  }

  //! Construct (RhoxD_E) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesRhoxDEOperator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesRhoxDT2B0WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxDT2B1WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxDT2Rho1WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxDT2B2WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxDT2Rho2WaveFunc(wvfParams));

    return wvf;
  }



  // Construct (Pion_2xD_T2) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesPion2xDT2Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesPion2xDT2Pi0WaveFunc(wvfParams));
    wvf.push_back(new MesPion2xDT2Pi1WaveFunc(wvfParams));
    wvf.push_back(new MesPion2xDT2A1WaveFunc(wvfParams));
    wvf.push_back(new MesPion2xDT2Pi2WaveFunc(wvfParams));
    wvf.push_back(new MesPion2xDT2A2WaveFunc(wvfParams));
    wvf.push_back(new MesPion2xDT2A3WaveFunc(wvfParams));

    return wvf;
  }

  // Construct (PionxD_T2) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesPionxDT2Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesPionxDT2Pi0WaveFunc(wvfParams));
    wvf.push_back(new MesPionxDT2A1WaveFunc(wvfParams));
    wvf.push_back(new MesPionxDT2Pi2WaveFunc(wvfParams));

    return wvf;
  }

  // Construct (A0xD_T2) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA0xDT2Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA0xDT2A0WaveFunc(wvfParams));
    wvf.push_back(new MesA0xDT2A2WaveFunc(wvfParams));
    wvf.push_back(new MesA0xDT2Pi1WaveFunc(wvfParams));

    return wvf;
  }
 
  // Construct (A0xB_T1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA0xBT1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA0xBT1B1WaveFunc(wvfParams));
    wvf.push_back(new MesA0xBT1Rho1WaveFunc(wvfParams));

    return wvf;
  }

  //! Construct (PionxB_T1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesPionxBT1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesPionxBT1Rho1WaveFunc(wvfParams));
    wvf.push_back(new MesPionxBT1B1WaveFunc(wvfParams));

    return wvf;
  }

  // Construct (A0_2xB_T1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA02xBT1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA02xBT1Pi0WaveFunc(wvfParams));
    wvf.push_back(new MesA02xBT1A1WaveFunc(wvfParams));
    wvf.push_back(new MesA02xBT1Pi1WaveFunc(wvfParams));
    wvf.push_back(new MesA02xBT1A2WaveFunc(wvfParams));
    wvf.push_back(new MesA02xBT1Pi2WaveFunc(wvfParams));

    return wvf;
  }

  // Construct (RhoxB_A1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesRhoxBA1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesRhoxBA1Pi0WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxBA1A1WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxBA1Pi2WaveFunc(wvfParams));

    return wvf;
  }

  // Construct (RhoxB_T1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesRhoxBT1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesRhoxBT1A0WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxBT1A1WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxBT1Pi1WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxBT1A2WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxBT1Pi2WaveFunc(wvfParams));

    return wvf;
  }


  //! Construct (RhoxB_T2) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesRhoxBT2Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesRhoxBT2A1WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxBT2Pi1WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxBT2A2WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxBT2Pi2WaveFunc(wvfParams));

    return wvf;
  }

  //! Construct (RhoxB_E) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesRhoxBEOperator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesRhoxBEA1WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxBEPi1WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxBEA2WaveFunc(wvfParams));
    wvf.push_back(new MesRhoxBEPi2WaveFunc(wvfParams));

    return wvf;
  }

  //! Construct (A1xB_A1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA1xBA1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA1xBA1Rho1WaveFunc(wvfParams));
    wvf.push_back(new MesA1xBA1B2WaveFunc(wvfParams));

    return wvf;
  }


  //! Construct (A1xB_T1) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA1xBT1Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA1xBT1Rho0WaveFunc(wvfParams));
    wvf.push_back(new MesA1xBT1Rho1WaveFunc(wvfParams));
    wvf.push_back(new MesA1xBT1B1WaveFunc(wvfParams));
    wvf.push_back(new MesA1xBT1Rho2WaveFunc(wvfParams));
    wvf.push_back(new MesA1xBT1B2WaveFunc(wvfParams));

    return wvf;
  }

  //! Construct (A1xB_T2) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA1xBT2Operator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA1xBT2Rho1WaveFunc(wvfParams));
    wvf.push_back(new MesA1xBT2B1WaveFunc(wvfParams));
    wvf.push_back(new MesA1xBT2Rho2WaveFunc(wvfParams));
    wvf.push_back(new MesA1xBT2B2WaveFunc(wvfParams));

    return wvf;
  }

  //! Construct (A1xB_E) operator
  std::list< Handle< WaveFunction<Complex> > > 
  MesA1xBEOperator::overlaps() const
  {
    std::list< Handle< WaveFunction<Complex> > > wvf;

    WaveFuncParams wvfParams;
    wvfParams.lattice = params.lattice;

    wvf.push_back(new MesA1xBERho1WaveFunc(wvfParams));
    wvf.push_back(new MesA1xBEB1WaveFunc(wvfParams));
    wvf.push_back(new MesA1xBERho2WaveFunc(wvfParams));
    wvf.push_back(new MesA1xBEB2WaveFunc(wvfParams));

    return wvf;
  }

  
  //! Meson sources
  namespace MesonOperatorEnv
  { 

    // Register all the possible deriv mesons
    bool registerAll(void) 
    {
      bool foo = true;

      //! Register all the factories
      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("Gamma0"),
								    mesGamma0Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("Gamma1"),
								    mesGamma1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("Gamma2"),
								    mesGamma2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("Gamma3"),
								    mesGamma3Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("Gamma4"),
								    mesGamma4Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("Gamma5"),
								    mesGamma5Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("Gamma6"),
								    mesGamma6Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("Gamma7"),
								    mesGamma7Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("Gamma8"),
								    mesGamma8Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("Gamma9"),
								    mesGamma9Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("Gamma10"),
								    mesGamma10Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("Gamma11"),
								    mesGamma11Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("Gamma12"),
								    mesGamma12Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("Gamma13"),
								    mesGamma13Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("Gamma14"),
								    mesGamma14Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("Gamma15"),
								    mesGamma15Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a0"),
								    mesA0A1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("pion"),
								    mesPionA1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("rho"),
								    mesRhoT1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a1"),
								    mesA1T1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("b1"),
								    mesB1T1Operator);


      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("pionxNABLA_T1"),
								    mesPionxNablaT1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a0xNABLA_T1"),
								    mesA0xNablaT1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("b0xNABLA_T1"),
								    mesA02xNablaT1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("rhoxNABLA_A1"),
								    mesRhoxNablaA1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("rhoxNABLA_T1"),
								    mesRhoxNablaT1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("rhoxNABLA_T2"),
								    mesRhoxNablaT2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("rhoxNABLA_E"),
								    mesRhoxNablaEOperator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a1xNABLA_A1"),
								    mesA1xNablaA1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a1xNABLA_T1"),
								    mesA1xNablaT2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a1xNABLA_T2"),
								    mesA1xNablaT2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a1xNABLA_E"),
								    mesA1xNablaEOperator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("b1xNABLA_A1"),
								    mesB1xNablaA1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("b1xNABLA_T1"),
								    mesB1xNablaT1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("b1xNABLA_T2"),
								    mesB1xNablaT2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("b1xNABLA_E"),
								    mesB1xNablaEOperator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("b0xD_T2"),
								    mesA02xDT2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a1xD_A2"),
								    mesA1xDA2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a1xD_E"),
								    mesA1xDEOperator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a1xD_T1"),
								    mesA1xDT1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a1xD_T2"),
								    mesA1xDT2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("b1xD_A2"),
								    mesB1xDA2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("b1xD_E"),
								    mesB1xDEOperator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("b1xD_T1"),
								    mesB1xDT1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("b1xD_T2"),
								    mesB1xDT2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("rhoxD_A2"),
								    mesRhoxDA2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("rhoxD_T1"),
								    mesRhoxDT1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("rhoxD_T2"),
								    mesRhoxDT2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("rhoxD_E"),
								    mesRhoxDEOperator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("pion_2xD_T2"),
								    mesPion2xDT2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("pionxD_T2"),
								    mesPionxDT2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a0xD_T2"),
								    mesA0xDT2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a0xB_T1"),
								    mesA0xBT1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("pionxB_T1"),
								    mesPionxBT1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("b0xB_T1"),
								    mesA02xBT1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("rhoxB_A1"),
								    mesRhoxBA1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("rhoxB_T1"),
								    mesRhoxBT1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("rhoxB_T2"),
								    mesRhoxBT2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("rhoxB_E"),
								    mesRhoxBEOperator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a1xB_A1"),
								    mesA1xBA1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a1xB_T1"),
								    mesA1xBT1Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a1xB_T2"),
								    mesA1xBT2Operator);

      foo &= FF::TheMesonOperatorFactory::Instance().registerObject(std::string("a1xB_E"),
								    mesA1xBEOperator);

      return foo;
    }

    const bool registered = registerAll();

  }  // end namespace MesonOperatorEnv

}  // end namespace FF



  
