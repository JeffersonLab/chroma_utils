// -*- C++ -*-
// $Id: formfac_operators.h,v 2.0 2008/12/05 04:43:36 edwards Exp $

/*! \file
 * \brief Operators and wavefunction overlaps
 */

#ifndef __operators_h__
#define __operators_h__

#include "adat/handle.h"
#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
// #include "formfac/formfac.h"
#include "wavefuncs.h"
#include <list>

namespace FF
{
  using ADAT::Handle;

  //----------------------------------------------------------------------------------
  //! Base class for Minkowski-space operators
  template<typename T>
  class Operator
  {
  public:
    //! The name of this interpolating field
    virtual std::string operatorName() const = 0;

    //! The allowed wavefunction overlaps 
    virtual std::list< Handle< WaveFunction<T> > > overlaps() const = 0;
  };


  //----------------------------------------------------------------------------------
  
  //! Params for derivative quark displacement
  struct OperatorParams
  {
    OperatorParams();
    OperatorParams(XMLReader& in, const std::string& path);
    void writeXML(XMLWriter& in, const std::string& path) const;

    std::string        operator_type;      /*!< Id of operator */

    LatticeParam       lattice;            /*!< Holds lattice size and aniso*/
  };



  //----------------------------------------------------------------------------------
  // Simple \Gamma(n) operators

  //! Construct (Gamma(0)) operator
  /*!
   * Operator is  Gamma(0)
   * The interpolator structure is
   * \f$\Gamma(0)\f$
   */
  class MesGamma0Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesGamma0Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (Gamma(1)) operator
  /*!
   * Operator is  Gamma(1)
   * The interpolator structure is
   * \f$\Gamma(1)\f$
   */
  class MesGamma1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesGamma1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return "gamma1";}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (Gamma(2)) operator
  /*!
   * Operator is  Gamma(2)
   * The interpolator structure is
   * \f$\Gamma(2)\f$
   */
  class MesGamma2Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesGamma2Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (Gamma(3)) operator
  /*!
   * Operator is  Gamma(3)
   * The interpolator structure is
   * \f$\Gamma(3)\f$
   */
  class MesGamma3Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesGamma3Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (Gamma(4)) operator
  /*!
   * Operator is  Gamma(4)
   * The interpolator structure is
   * \f$\Gamma(4)\f$
   */
  class MesGamma4Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesGamma4Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (Gamma(5)) operator
  /*!
   * Operator is  Gamma(5)
   * The interpolator structure is
   * \f$\Gamma(5)\f$
   */
  class MesGamma5Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesGamma5Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (Gamma(6)) operator
  /*!
   * Operator is  Gamma(6)
   * The interpolator structure is
   * \f$\Gamma(6)\f$
   */
  class MesGamma6Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesGamma6Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (Gamma(7)) operator
  /*!
   * Operator is  Gamma(7)
   * The interpolator structure is
   * \f$\Gamma(7)\f$
   */
  class MesGamma7Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesGamma7Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (Gamma(8)) operator
  /*!
   * Operator is  Gamma(8)
   * The interpolator structure is
   * \f$\Gamma(8)\f$
   */
  class MesGamma8Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesGamma8Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (Gamma(9)) operator
  /*!
   * Operator is  Gamma(9)
   * The interpolator structure is
   * \f$\Gamma(9)\f$
   */
  class MesGamma9Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesGamma9Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (Gamma(10)) operator
  /*!
   * Operator is  Gamma(10)
   * The interpolator structure is
   * \f$\Gamma(10)\f$
   */
  class MesGamma10Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesGamma10Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (Gamma(11)) operator
  /*!
   * Operator is  Gamma(11)
   * The interpolator structure is
   * \f$\Gamma(11)\f$
   */
  class MesGamma11Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesGamma11Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (Gamma(12)) operator
  /*!
   * Operator is  Gamma(12)
   * The interpolator structure is
   * \f$\Gamma(12)\f$
   */
  class MesGamma12Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesGamma12Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (Gamma(13)) operator
  /*!
   * Operator is  Gamma(13)
   * The interpolator structure is
   * \f$\Gamma(13)\f$
   */
  class MesGamma13Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesGamma13Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (Gamma(14)) operator
  /*!
   * Operator is  Gamma(14)
   * The interpolator structure is
   * \f$\Gamma(14)\f$
   */
  class MesGamma14Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesGamma14Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (Gamma(15)) operator
  /*!
   * Operator is  Gamma(15)
   * The interpolator structure is
   * \f$\Gamma(15)\f$
   */
  class MesGamma15Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesGamma15Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (a_0) operator
  /*!
   * Operator is  a_0
   * The interpolator structure is
   * \f$1\f$
   */
  class MesA0A1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA0A1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (pi_0) operator
  /*!
   * Operator is  pi_0
   * The interpolator structure is
   * \f$i\gamma^5\f$
   */
  class MesPionA1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesPionA1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (rho_1) operator
  /*!
   * Operator is  rho_1
   * The interpolator structure is
   * \f$\gamma^i\f$
   */
  class MesRhoT1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesRhoT1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (a_1) operator
  /*!
   * Operator is  a_1
   * The interpolator structure is
   * \f$\gamma^5\gamma^i\f$
   */
  class MesA1T1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA1T1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (b_1) operator
  /*!
   * Operator is  b_1
   * The interpolator structure is
   * \f$\sigma^{\mu\nu}\f$
   */
  class MesB1T1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesB1T1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (pionxNABLA_T1) operator
  /*!
   * Operator is  pion x Nabla_T1
   * The interpolator structure is
   * \f$\gamma_5\nabla_i\f$
   */
  class MesPionxNablaT1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesPionxNablaT1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (a0xaNABLA_T1) operator
  /*!
   * Operator is  a0 x nabla_T1
   * The interpolator is   
   * \f$\nabla_i\f$
   */
  class MesA0xNablaT1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA0xNablaT1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (b0xNABLA_T1) operator
  /*!
   * Operator is  b0 x nabla_T1
   * The interpolator is   
   * \f$\gamma_4 \nabla_i\f$
   */
  class MesA02xNablaT1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA02xNablaT1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (RhoxNABLA_A1) operator
  /*!
   * Operator is  rho x nabla_A1
   * The interpolator is   
   * \f$\gamma_i\nabla_i\f$
   */
  class MesRhoxNablaA1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesRhoxNablaA1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (rhoxNABLA_T1) operator
  /*!
   * Operator is  rho x nabla_T1
   * The interpolator is   
   * \f$\epsilon_{ijk}\gamma_j \nabla_k\f$  
   */
  class MesRhoxNablaT1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesRhoxNablaT1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (RhoxNABLA_T2) operator
  /*!
   * Operator is  rho x nabla_T2
   * The interpolator is   
   * \f$s_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxNablaT2Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesRhoxNablaT2Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

  //! Construct (RhoxNABLA_E) operator
  /*!
   * Operator is  rho x nabla_E
   * The interpolator is   
   * \f$S_{ajk}\gamma_j D_k\f$  
   */
  class MesRhoxNablaEOperator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesRhoxNablaEOperator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (a1xNABLA_A1) operator
  /*!
   * Operator is  a1 x D_A1
   * The interpolator is   
   * \f$\gamma_5\gamma_i \nabla_i\f$  
   */
  class MesA1xNablaA1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA1xNablaA1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

 //! Construct (a1xNABLA_T1) operator
  /*!
   * Operator is  a1 x nabla_T1
   * The interpolator is   
   * \f$\gamma_5 \epsilon_{ijk}\gamma_j \nabla_k\f$  
   */
  class MesA1xNablaT1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA1xNablaT1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

  //! Construct (a1xNABLA_T2) operator
  /*!
   * Operator is  a1 x nabla_T2
   * The interpolator is   
   * \f$\gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
   */
  class MesA1xNablaT2Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA1xNablaT2Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (A1xNABLA_E) operator
  /*!
   * Operator is  a1 x nabla_E
   * The interpolator is   
   * \f$\gamma_5 S_{\alpha jk}\gamma_j \nabla_k\f$  
   */
  class MesA1xNablaEOperator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA1xNablaEOperator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

  //! Construct (b1xNABLA_A1) operator
  /*!
   * Operator is  b1 x nabla_A1
   * The interpolator is   
   * \f$\gamma_4\gamma_5\gamma_k \nabla_k\f$  
   */
  class MesB1xNablaA1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaA1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

  //! Construct (b1xNABLA_T1) operator
  /*!
   * Operator is  b1 x nabla_T1
   * The interpolator is   
   * \f$\gamma_4\gamma_5\epsilon_{ijk}\gamma_j \nabla_k\f$  
   */
  class MesB1xNablaT1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaT1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

  //! Construct (b1xNABLA_T2) operator
  /*!
   * Operator is  b1 x nabla_T2
   * The interpolator is   
   * \f$\gamma_4\gamma_5s_{ijk}\gamma_j \nabla_k\f$  
   */
  class MesB1xNablaT2Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaT2Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

  //! Construct (b1xNABLA_E) operator
  /*!
   * Operator is  b1 x nabla_TE
   * The interpolator is   
   * \f$\gamma_4\gamma_5 S_{ijk}\gamma_j \nabla_k\f$  
   */
  class MesB1xNablaEOperator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaEOperator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (b0xD_T2) operator
  /*!
   * Operator is  b0 x D_T2
   * The interpolator is   
   * \f$\gamma_4 D_i\f$  
   */
  class MesA02xDT2Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA02xDT2Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (a1xD_A2) operator
  /*!
   * Operator is  a1 x D_A2
   * The interpolator is   
   * \f$\gamma_5\gamma_i D_i\f$  
   */
  class MesA1xDA2Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA1xDA2Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (a1xD_E) operator
  /*!
   * Operator is  a1 x D_E
   * The interpolator is   
   * \f$\gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
   */
  class MesA1xDEOperator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA1xDEOperator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (a1xD_T1) operator
  /*!
   * Operator is  a1 x D_T1
   * The interpolator is   
   * \f$\gamma_5 s_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDT1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA1xDT1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (a1xD_T2) operator
  /*!
   * Operator is  a1 x D_T2
   * The interpolator is   
   * \f$\gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDT2Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA1xDT2Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (b1xD_A2) operator
  /*!
   * Operator is  b1 x D_A2
   * The interpolator is   
   * \f$\gamma_4\gamma_5 \gamma_i D_i\f$  
   */
  class MesB1xDA2Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesB1xDA2Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

  //! Construct (B1xD_E) operator
  /*!
   * Operator is  b1 x D_E
   * The interpolator is   
   * \f$\gamma_4\gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
   */
  class MesB1xDEOperator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesB1xDEOperator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (B1xD_T1) operator
  /*!
   * Operator is  b1 x D_T1
   * The interpolator is   
   * \f$\gamma_4\gamma_5 s_{ijk}\gamma_j D_k\f$  
   */
  class MesB1xDT1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesB1xDT1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (B1xD_T2) operator
  /*!
   * Operator is  b1 x D_T2
   * The interpolator is   
   * \f$\gamma_4\gamma_5 \epsilon_{ijk}\gamma_j D_k\f$  
   */
  class MesB1xDT2Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesB1xDT2Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (RhoxD_A2) operator
  /*!
   * Operator is  rho x D_A2
   * The interpolator is   
   * \f$\gamma_i D_i\f$  
   */
  class MesRhoxDA2Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDA2Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (RhoxD_T1) operator
  /*!
   * Operator is  rho x D_T1
   * The interpolator is   
   * \f$s_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDT1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDT1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (RhoxD_T2) operator
  /*!
   * Operator is  rho x D_T2
   * The interpolator is   
   * \f$\epsilon_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDT2Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDT2Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

  //! Construct (RhoxD_E) operator
  /*!
   * Operator is  rho x D_E
   * The interpolator is   
   * \f$S_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDEOperator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDEOperator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (pion_2xD_T2) operator
  /*!
   * Operator is  pion_2 x D_T2
   * The interpolator is   
   * \f$\gamma_4\gamma_5 D_i\f$  
   */
  class MesPion2xDT2Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesPion2xDT2Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

  //! Construct (pionxD_T2) operator
  /*!
   * Operator is  pion x D_T2
   * The interpolator is   
   * \f$\gamma_5 D_i\f$  
   */
  class MesPionxDT2Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesPionxDT2Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

  //! Construct (a0xD_T2) operator
  /*!
   * Operator is  a0 x D_T2
   * The interpolator is   
   * \f$D_i\f$  
   */
  class MesA0xDT2Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA0xDT2Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

  //! Construct (a0xB_T1) operator
  /*!
   * Operator is  a0 x B_T1
   * The interpolator is   
   * \f$B_i\f$  
   */
  class MesA0xBT1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA0xBT1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (pionxB_T1) operator
  /*!
   * Operator is  pion x B_T1
   * The interpolator is   
   * \f$\gamma_5 B_i\f$  
   */
  class MesPionxBT1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesPionxBT1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

  //! Construct (b0xB_T1) operator
  /*!
   * Operator is  b0 x B_T1
   * The interpolator is   
   * \f$\gamma_4  B_i\f$  
   */
  class MesA02xBT1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA02xBT1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

  //! Construct (RhoxB_TA) operator
  /*!
   * Operator is rho x B_A1
   * The interpolator is   
   * \f$\gamma_k B_k\f$  
   */
  class MesRhoxBA1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBA1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (RhoxB_T1) operator
  /*!
   * Operator is rho x B_T1
   * The interpolator is   
   * \f$\epsilon_{ijk}\gamma_j B_k\f$  
   */
  class MesRhoxBT1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBT1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (RhoxB_T2) operator
  /*!
   * Operator is  rho x B_T2
   * The interpolator is   
   * \f$s_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxBT2Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBT2Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

  //! Construct (RhoxB_E) operator
  /*!
   * Operator is  rho x B_E
   * The interpolator is   
   * \f$S_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxBEOperator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBEOperator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

  //! Construct (a1xB_A1) operator
  /*!
   * Operator is  a1 x B_A1
   * The interpolator is   
   * \f$\gamma_5 \gamma_i B_i\f$  
   */
  class MesA1xBA1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA1xBA1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (a1xB_T1) operator
  /*!
   * Operator is  a1 x B_T1
   * The interpolator is   
   * \f$\gamma_5 \epsilon_{ijk}\gamma_j B_k\f$  
   */
  class MesA1xBT1Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA1xBT1Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Construct (A1xB_T2) operator
  /*!
   * Operator is  a1 x B_T2
   * The interpolator is   
   * \f$\gamma_5 s_{ijk}\gamma_j B_k\f$  
   */
  class MesA1xBT2Operator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA1xBT2Operator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };

 //! Construct (A1xB_E) operator
  /*!
   * Operator is  a1 x B_E
   * The interpolator is   
   * \f$\gamma_5 S_{ijk}\gamma_j B_k\f$  
   */
  class MesA1xBEOperator : public Operator<Complex>
  {
  public:
    //! Full constructor
    MesA1xBEOperator(const OperatorParams& p) : params(p) {}

    //! The name of this interpolating field
    std::string operatorName() const {return params.operator_type;}

    //! The allowed wavefunction overlaps 
    std::list< Handle< WaveFunction<Complex> > > overlaps() const;

  private:
    OperatorParams  params;   /*!< operator params */
  };


  //! Reader
  void read(XMLReader& xml, const std::string& path, OperatorParams& param);

  //! Writer
  void write(XMLWriter& xml, const std::string& path, const OperatorParams& param);

} // namespace FF

#endif
