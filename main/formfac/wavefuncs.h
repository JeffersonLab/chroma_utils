// -*- C++ -*-
// $Id: wavefuncs.h,v 2.0 2008/12/05 04:43:49 edwards Exp $

/*! \file
 * \brief Wavefunction overlaps
 */

#ifndef __wavefuncs_h__
#define __wavefuncs_h__

#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "formfac/formfac_ensemble.h"

namespace FF
{

  Array<Complex> minkPolVec(const Real& mass, const ArrayInt& p, int r, LatticeParam lattice);

  //  Array<Complex> minkPolVecDUM(const Real& mass, const ArrayInt& p, int r, LatticeParam lattice);

  Array< Array<Complex> > minkPolTens(const Real& mass, const ArrayInt& p, int r, LatticeParam lattice);

  //  Array< Array<Complex> > minkPolTensDUM(const Real& mass, const ArrayInt& p, int r, LatticeParam lattice);

  Array<Complex> test_wvf0index(const Real& mass, const ArrayInt& p, LatticeParam lattice);
  Array<Complex> test_wvf1index(const Real& mass, const ArrayInt& p, int r, LatticeParam lattice, int mu);
  Array<Complex> test_wvf2index(const Real& mass, const ArrayInt& p, int r, LatticeParam lattice, int mu, int nu);
  Real test_antisym();
  Real test_antisym3();

  Real cleb(double j1, double m1, double j2, double m2, double j3, double m3);

  //! Make a Minkowski four vector momentum
  /*! \return \f$p^{\mu}\f$*/
  Array<Real> make4Vec(const Real& mass, const ArrayInt& p, LatticeParam lattice);


  //----------------------------------------------------------------------------------
  //! Base class for Minkowski-space wavefunction overlaps
  template<typename T>
  class WaveFunction
  {
  public:
    //! The name of this particle
    virtual std::string particleName() const = 0;

    //! The name of this wavefunction
//    virtual std::string waveFuncName() const = 0;

    //! The number of polarizations
    virtual int numPolar() const = 0;

    //! The number of directions
    virtual int numDir() const = 0;

    //! The overlap factor
    /*! 
     *  \param p    Minkowski-space 4-vector
     *  \param dir  Free direction-like index (1-based)
     *  \param r    Polarization index (0-based)
     */
    virtual Array<T> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const = 0;
  };


  //----------------------------------------------------------------------------------
  //! Params for derivative quark displacement
  struct WaveFuncParams
  {
    WaveFuncParams();
    WaveFuncParams(XMLReader& in, const std::string& path);
    void writeXML(XMLWriter& in, const std::string& path) const;

    std::string        wavefunc_type;      /*!< Id of wavefunc */

    LatticeParam       lattice;            /*!< Holds lattice size and aniso*/
  };


  //! Reader
  void read(XMLReader& xml, const std::string& path, WaveFuncParams& param);

  //! Writer
  void write(XMLWriter& xml, const std::string& path, const WaveFuncParams& param);



  //----------------------------------------------------------------------------------
  //! Construct (Gamma(0)) wavefunction overlap onto an  a_0
  /*!
   * WaveFunction is  Gamma(0)  on an  a_0  state
   * The interpolator structure is
   * \f$\Gamma(0)\f$
   */
  class MesGamma0A0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma0A0WaveFunc(const WaveFuncParams& p) : params(p) {      }

    //! The name of this particle
    std::string particleName() const {return "a0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The number of kinematic factors
    //int numTerms() const {return len;} 

    //! The overlap factor - now an array
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
    //int len;
  };


  //----------------------------------------------------------------------------------
  //! Construct (Gamma(1)) wavefunction overlap onto an  b_0
  /*!
   * WaveFunction is  Gamma(1)  on an  b_0  state
   * The interpolator structure is
   * \f$\Gamma(1)\f$
   */
  class MesGamma1B0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma1B0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (Gamma(1)) wavefunction overlap onto an  rho_1
  /*!
   * WaveFunction is  Gamma(1) on an  rho_1  state
   * The interpolator structure is
   * \f$\Gamma(1)\f$
   */
  class MesGamma1Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma1Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (Gamma(2)) wavefunction overlap onto an  b_0
  /*!
   * WaveFunction is  Gamma(1)  on an  b_0  state
   * The interpolator structure is
   * \f$\Gamma(1)\f$
   */
  class MesGamma2B0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma2B0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (Gamma(2)) wavefunction overlap onto an  rho_1
  /*!
   * WaveFunction is  Gamma(2) on an  rho_1  state
   * The interpolator structure is
   * \f$\Gamma(2)\f$
   */
  class MesGamma2Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma2Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (Gamma(3)) wavefunction overlap onto an  b_1
  /*!
   * WaveFunction is  Gamma(3)  on an  b_1  state
   * The interpolator structure is
   * \f$\Gamma(3)\f$
   */
  class MesGamma3B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma3B1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (Gamma(3)) wavefunction overlap onto an  rho_1
  /*!
   * WaveFunction is  Gamma(3) on an  rho_1  state
   * The interpolator structure is
   * \f$\Gamma(3)\f$
   */
  class MesGamma3Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma3Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (Gamma(4)) wavefunction overlap onto an  b_0
  /*!
   * WaveFunction is  Gamma(4)  on an  b_0  state
   * The interpolator structure is
   * \f$\Gamma(4)\f$
   */
  class MesGamma4B0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma4B0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (Gamma(4)) wavefunction overlap onto an  rho_1
  /*!
   * WaveFunction is  Gamma(4) on an  rho_1  state
   * The interpolator structure is
   * \f$\Gamma(4)\f$
   */
  class MesGamma4Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma4Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (Gamma(5)) wavefunction overlap onto an  b_1
  /*!
   * WaveFunction is  Gamma(5)  on an  b_1  state
   * The interpolator structure is
   * \f$\Gamma(5)\f$
   */
  class MesGamma5B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma5B1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (Gamma(5)) wavefunction overlap onto an  rho_1
  /*!
   * WaveFunction is  Gamma(5) on an  rho_1  state
   * The interpolator structure is
   * \f$\Gamma(5)\f$
   */
  class MesGamma5Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma5Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (Gamma(6)) wavefunction overlap onto an  b_1
  /*!
   * WaveFunction is  Gamma(6)  on an  b_1  state
   * The interpolator structure is
   * \f$\Gamma(6)\f$
   */
  class MesGamma6B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma6B1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (Gamma(6)) wavefunction overlap onto an  rho_1
  /*!
   * WaveFunction is  Gamma(6) on an  rho_1  state
   * The interpolator structure is
   * \f$\Gamma(6)\f$
   */
  class MesGamma6Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma6Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (Gamma(7)) wavefunction overlap onto an  pi_0
  /*!
   * WaveFunction is  Gamma(7)  on an  pi_0  state
   * The interpolator structure is
   * \f$\Gamma(7)\f$
   */
  class MesGamma7Pi0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma7Pi0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (Gamma(7)) wavefunction overlap onto an  a_1
  /*!
   * WaveFunction is  Gamma(7) on an  a_1  state
   * The interpolator structure is
   * \f$\Gamma(7)\f$
   */
  class MesGamma7A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma7A1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (Gamma(8)) wavefunction overlap onto an  b_0
  /*!
   * WaveFunction is  Gamma(8)  on an  b_0  state
   * The interpolator structure is
   * \f$\Gamma(8)\f$
   */
  class MesGamma8B0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma8B0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (Gamma(8)) wavefunction overlap onto an  rho_1
  /*!
   * WaveFunction is  Gamma(8) on an  rho_1  state
   * The interpolator structure is
   * \f$\Gamma(8)\f$
   */
  class MesGamma8Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma8Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (Gamma(9)) wavefunction overlap onto an  b_1
  /*!
   * WaveFunction is  Gamma(9)  on an  b_1  state
   * The interpolator structure is
   * \f$\Gamma(9)\f$
   */
  class MesGamma9B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma9B1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (Gamma(9)) wavefunction overlap onto an  rho_1
  /*!
   * WaveFunction is  Gamma(9) on an  rho_1  state
   * The interpolator structure is
   * \f$\Gamma(9)\f$
   */
  class MesGamma9Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma9Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (Gamma(10)) wavefunction overlap onto an  b_1
  /*!
   * WaveFunction is  Gamma(10)  on an  b_1  state
   * The interpolator structure is
   * \f$\Gamma(10)\f$
   */
  class MesGamma10B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma10B1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (Gamma(10)) wavefunction overlap onto an  rho_1
  /*!
   * WaveFunction is  Gamma(10) on an  rho_1  state
   * The interpolator structure is
   * \f$\Gamma(10)\f$
   */
  class MesGamma10Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma10Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (Gamma(11)) wavefunction overlap onto an  pi_0
  /*!
   * WaveFunction is  Gamma(11)  on an  pi_0  state
   * The interpolator structure is
   * \f$\Gamma(11)\f$
   */
  class MesGamma11Pi0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma11Pi0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (Gamma(11)) wavefunction overlap onto an  a_1
  /*!
   * WaveFunction is  Gamma(11) on an  a_1  state
   * The interpolator structure is
   * \f$\Gamma(11)\f$
   */
  class MesGamma11A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma11A1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (Gamma(12)) wavefunction overlap onto an  b_1
  /*!
   * WaveFunction is  Gamma(12)  on an  b_1  state
   * The interpolator structure is
   * \f$\Gamma(12)\f$
   */
  class MesGamma12B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma12B1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (Gamma(12)) wavefunction overlap onto an  rho_1
  /*!
   * WaveFunction is  Gamma(12) on an  rho_1  state
   * The interpolator structure is
   * \f$\Gamma(12)\f$
   */
  class MesGamma12Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma12Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (Gamma(13)) wavefunction overlap onto an  pi_0
  /*!
   * WaveFunction is  Gamma(13)  on an  pi_0  state
   * The interpolator structure is
   * \f$\Gamma(13)\f$
   */
  class MesGamma13Pi0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma13Pi0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (Gamma(13)) wavefunction overlap onto an  a_1
  /*!
   * WaveFunction is  Gamma(13) on an  a_1  state
   * The interpolator structure is
   * \f$\Gamma(13)\f$
   */
  class MesGamma13A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma13A1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (Gamma(14)) wavefunction overlap onto an  pi_0
  /*!
   * WaveFunction is  Gamma(14)  on an  pi_0  state
   * The interpolator structure is
   * \f$\Gamma(14)\f$
   */
  class MesGamma14Pi0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma14Pi0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (Gamma(14)) wavefunction overlap onto an  a_1
  /*!
   * WaveFunction is  Gamma(14) on an  a_1  state
   * The interpolator structure is
   * \f$\Gamma(14)\f$
   */
  class MesGamma14A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma14A1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (Gamma(15)) wavefunction overlap onto an  pi_0
  /*!
   * WaveFunction is  Gamma(15)  on an  pi_0  state
   * The interpolator structure is
   * \f$\Gamma(15)\f$
   */
  class MesGamma15Pi0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesGamma15Pi0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (a_0) wavefunction overlap onto an  a_0
  /*!
   * WaveFunction is  1  on an  a_0  state
   * The interpolator structure is
   * \f$1\f$
   */
  class MesA0A1A0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA0A1A0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (pion) wavefunction overlap onto an  pi
  /*!
   * WaveFunction is  gamma^5  on an  pi  state
   * The interpolator structure is
   * \f$i\gamma^5\f$
   */
  class MesPionA1Pi0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesPionA1Pi0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (rho) wavefunction overlap onto an  b_0
  /*!
   * WaveFunction is  gamma^i  on an  b_0  state
   * The interpolator structure is
   * \f$gamma^i\f$
   */
  class MesRhoT1B0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoT1B0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (rho) wavefunction overlap onto an  rho_1
  /*!
   * WaveFunction is  gamma^i  on an  rho_1  state
   * The interpolator structure is
   * \f$gamma^i\f$
   */
  class MesRhoT1Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoT1Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (a_1) wavefunction overlap onto an  pi_0
  /*!
   * WaveFunction is  gamma_5 gamma_i  on an  pi_0  state
   * The interpolator structure is
   * \f$gamma^5\gamma^i\f$
   */
  class MesA1T1Pi0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1T1Pi0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (a_1) wavefunction overlap onto an  a_1
  /*!
   * WaveFunction is  gamma_5 gamma_i  on an  a_1  state
   * The interpolator structure is
   * \f$gamma^5\gamma^i\f$
   */
  class MesA1T1A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1T1A1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (b_1) wavefunction overlap onto an  b_1
  /*!
   * WaveFunction is  sigma^{mu,nu}  on an  b_1  state
   * The interpolator structure is
   * \f$\sigma^{\mu\nu}\f$
   */
  class MesB1T1B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1T1B1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (b_1) wavefunction overlap onto an  rho_1
  /*!
   * WaveFunction is  sigma^{mu,nu}  on an  rho_1  state
   * The interpolator structure is
   * \f$\sigma^{\mu\nu}\f$
   */
  class MesB1T1Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1T1Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (PionxNabla_T1) wavefunction overlap onto an  rho_0
  /*!
   * WaveFunction is  Pion x Nabla_T1 on an  rho_0  state
   * The interpolator structure is
   * \f$\gamma_5\nabla_i\f$
   */
  class MesPionxNablaT1Rho0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesPionxNablaT1Rho0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (PionxNabla_T1) wavefunction overlap onto a  b_1
  /*!
   * WaveFunction is  Pion x Nabla_T1 on a  b_1  state
   * The sink interpolator structure is
   * \f$\gamma_5\nabla_i\f$
   */
  class MesPionxNablaT1B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesPionxNablaT1B1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (A0xNabla_T1) wavefunction overlap onto a  b_0
  /*!
   * WaveFunction is  a0 x nabla_T1 on a b_0 state
   * The interpolator is   
   * \f$\nabla_i\f$
   */
  class MesA0xNablaT1B0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA0xNablaT1B0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (A0xNabla_T1) wavefunction overlap onto a  rho_1
  /*!
   * WaveFunction is  a0 x nabla_T1 on a rho_1 state
   * The interpolator is   
   * \f$\nabla_i\f$
   */
  class MesA0xNablaT1Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA0xNablaT1Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (A0_2xNabla_T1) wavefunction overlap onto a  a_0
  /*!
   * WaveFunction is  a0_2 x nabla_T1 on a a_0 state
   * The interpolator is   
   * \f$\gamma_4 \nabla_i\f$
   */
  class MesA02xNablaT1A0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA02xNablaT1A0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (A0_2xNabla_T1) wavefunction overlap onto a  a_1
  /*!
   * WaveFunction is  a0_2 x nabla_T1 on a a_1 state
   * The interpolator is   
   * \f$\gamma_4 \nabla_i\f$
   */
  class MesA02xNablaT1A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA02xNablaT1A1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (A0_2xNabla_T1) wavefunction overlap onto a  a_2
  /*!
   * WaveFunction is  a0_2 x nabla_T1 on a a_2 state
   * The interpolator is   
   * \f$\gamma_4 \nabla_i\f$
   */
  class MesA02xNablaT1A2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA02xNablaT1A2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (A0_2xNabla_T1) wavefunction overlap onto a  symmetric pi_1
  /*!
   * WaveFunction is  a0_2 x nabla_T1 on a pi_1 state
   * The interpolator is   
   * \f$\gamma_4 \nabla_i\f$
   */
  class MesA02xNablaT1Pi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA02xNablaT1Pi1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (RhoxNabla_A1) wavefunction onto an  a_0  
  /*!
   * WaveFunction is  rho x nabla_A1 on an  a_0 state
   * The interpolator is   
   * \f$\gamma_i\nabla_i\f$
   */
  class MesRhoxNablaA1A0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxNablaA1A0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //! Construct (RhoxNabla_A1) wavefunction onto an  a_2
  /*!
   * WaveFunction is  rho x nabla_A1 on an  a_2 state
   * The interpolator is   
   * \f$\gamma_i\nabla_i\f$
   */
  class MesRhoxNablaA1A2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxNablaA1A2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (RhoxNabla_A1) wavefunction onto an  pi_1
  /*!
   * WaveFunction is  rho x nabla_A1 on an  pi_1 state
   * The interpolator is   
   * \f$\gamma_i\nabla_i\f$
   */
  class MesRhoxNablaA1Pi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxNablaA1Pi1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (RhoxNabla_T1) wavefunction onto an a_1
  /*!
   * WaveFunction is  rho x nabla_T1 onto an a_1
   * The interpolator is   
   * \f$\epsilon_{ijk}\gamma_j \nabla_k\f$  
   */
  class MesRhoxNablaT1A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxNablaT1A1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //----------------------------------------------------------------------------------
  //! Construct (RhoxNabla_T1) wavefunction onto an \pi_1
  /*!
   * WaveFunction is  rho x nabla_T1 onto an \pi_1
   * The interpolator is   
   * \f$\epsilon_{ijk}\gamma_j \nabla_k\f$  
   */
  class MesRhoxNablaT1Pi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxNablaT1Pi1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //----------------------------------------------------------------------------------
  //! Construct (RhoxNabla_T2) wavefunction onto an a_0
  /*!
   * WaveFunction is  rho x nabla_T2 onto an a_0
   * The interpolator is   
   * \f$s_{ijk}\gamma_j \nabla_k\f$  
   */
  class MesRhoxNablaT2A0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxNablaT2A0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (RhoxNabla_T2) wavefunction onto an a_2
  /*!
   * WaveFunction is  rho x nabla_T2 onto an a_2
   * The interpolator is   
   * \f$s_{ijk}\gamma_j \nabla_k\f$  
   */
  class MesRhoxNablaT2A2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxNablaT2A2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (RhoxNabla_T2) wavefunction onto a \pi_1
  /*!
   * WaveFunction is  rho x nabla_T2 onto an \pi_1 
   * The interpolator is   
   * \f$s_{ijk}\gamma_j \nabla_k\f$  
   */
  class MesRhoxNablaT2Pi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxNablaT2Pi1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //----------------------------------------------------------------------------------
  //! Construct (RhoxNabla_E) wavefunction onto an a_0
  /*!
   * WaveFunction is  rho x nabla_E onto an a_0
   * The interpolator is   
   * \f$S_{\alpha jk}\gamma_j \nabla_k\f$  
   */
  class MesRhoxNablaEA0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxNablaEA0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //----------------------------------------------------------------------------------
  //! Construct (RhoxNabla_E) wavefunction onto an a_2
  /*!
   * WaveFunction is  rho x nabla_E onto an a_0 - pp-type
   * The interpolator is   
   * \f$S_{\alpha jk}\gamma_j \nabla_k\f$  
   */
  class MesRhoxNablaEA2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxNablaEA2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //----------------------------------------------------------------------------------
  //! Construct (RhoxNabla_E) wavefunction onto an \pi_1
  /*!
   * WaveFunction is  rho x nabla_E onto an \pi_1 
   * The interpolator is   
   * \f$S_{\alpha jk}\gamma_j \nabla_k\f$  
   */
  class MesRhoxNablaEPi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxNablaEPi1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //----------------------------------------------------------------------------------
  //! Construct (A1xNabla_A1) wavefunction onto a \rho_0
  /*!
   * WaveFunction is  a1 x D_A1 onto a \rho_0
   * The interpolator is   
   * \f$\gamma_5\gamma_i \nabla_i\f$  
   */
  class MesA1xNablaA1Rho0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xNablaA1Rho0WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (A1xNabla_A1) wavefunction onto a \rho_2
  /*!
   * WaveFunction is  a1 x D_A1 onto a \rho_2
   * The interpolator is   
   * \f$\gamma_5\gamma_i \nabla_i\f$  
   */
  class MesA1xNablaA1Rho2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xNablaA1Rho2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho0";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (A1xNabla_A1) wavefunction onto a b_1
  /*!
   * WaveFunction is  a1 x D_A1 onto a b_1
   * The interpolator is   
   * \f$\gamma_5\gamma_i \nabla_i\f$  
   */
  class MesA1xNablaA1B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xNablaA1B1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


 //----------------------------------------------------------------------------------
  //! Construct (A1xNabla_T1) wavefunction onto a \rho_1
  /*!
   * WaveFunction is  a1 x NABLA_T1 onto a \rho_1
   * The interpolator is   
   * \f$\epsilon_{ijk} \gamma_5\gamma_j \nabla_k\f$  
   */
  class MesA1xNablaT1Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xNablaT1Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (A1xNabla_T1) wavefunction onto a b_1
  /*!
   * WaveFunction is  a1 x NABLA_T1 onto a b_1 
   * The interpolator is   
   * \f$\epsilon_{ijk} \gamma_5\gamma_j \nabla_k\f$  
   */
  class MesA1xNablaT1B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xNablaT1B1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (A1xNabla_T2) wavefunction onto a \rho_0
  /*!
   * WaveFunction is  a1 x nabla_T2 - onto a \rho_0 
   * The interpolator is   
   * \f$\gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
   */
  class MesA1xNablaT2Rho0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xNablaT2Rho0WaveFunc(const WaveFuncParams& p) : params(p) {}
    
    //! The name of this particle
    std::string particleName() const {return "rho0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (A1xNabla_T2) wavefunction onto a \rho_2
  /*!
   * WaveFunction is  a1 x nabla_T2 - onto a \rho_2
   * The interpolator is   
   * \f$\gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
   */
  class MesA1xNablaT2Rho2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xNablaT2Rho2WaveFunc(const WaveFuncParams& p) : params(p) {}
    
    //! The name of this particle
    std::string particleName() const {return "rho2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (A1xNabla_T2) wavefunction onto a b_1
  /*!
   * WaveFunction is  a1 x nabla_T2 - onto a b_1 
   * The interpolator is   
   * \f$\gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
   */
  class MesA1xNablaT2B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xNablaT2B1WaveFunc(const WaveFuncParams& p) : params(p) {}
    
    //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (A1xNabla_E) wavefunction onto a \rho_0
  /*!
   * WaveFunction is  a1 x nabla_E onto a \rho_0 
   * The interpolator is   
   * \f$\gamma_5 S_{\alpha jk}\gamma_j \nabla_k\f$  
   */
  class MesA1xNablaERho0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xNablaERho0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (A1xNabla_E) wavefunction onto a \rho_2
  /*!
   * WaveFunction is  a1 x nabla_E onto a \rho_2
   * The interpolator is   
   * \f$\gamma_5 S_{\alpha jk}\gamma_j \nabla_k\f$  
   */
  class MesA1xNablaERho2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xNablaERho2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (A1xNabla_E) wavefunction onto a b_1
  /*!
   * WaveFunction is  a1 x nabla_E onto a b_1 
   * The interpolator is   
   * \f$\gamma_5 S_{\alpha jk}\gamma_j \nabla_k\f$  
   */
  class MesA1xNablaEB1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xNablaEB1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (B1xNabla_A1) wavefunction onto a \pi_0
  /*!
   * WaveFunction is  b1 x nabla_A1 onto a \pi_0 
   * The interpolator is   
   * \f$\gamma_4 \gamma_5 \gamma_j \nabla_j\f$  
   */
  class MesB1xNablaA1Pi0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaA1Pi0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //----------------------------------------------------------------------------------
  //! Construct (B1xNabla_A1) wavefunction onto a a1
  /*!
   * WaveFunction is  b1 x nabla_A1 onto a a1
   * The interpolator is   
   * \f$\gamma_4 \gamma_5 \gamma_j \nabla_j\f$  
   */
  class MesB1xNablaA1A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaA1A1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //----------------------------------------------------------------------------------
  //! Construct (B1xNabla_A1) wavefunction onto a \pi_2
  /*!
   * WaveFunction is  b1 x nabla_A1 onto a \pi_2
   * The interpolator is   
   * \f$\gamma_4 \gamma_5 \gamma_j \nabla_j\f$  
   */
  class MesB1xNablaA1Pi2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaA1Pi2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //----------------------------------------------------------------------------------
  //! Construct (B1xNabla_T1) wavefunction onto a a0
  /*!
   * WaveFunction is  b1 x nabla_T1 onto a a0
   * The interpolator is   
   * \f$\gamma_4 \gamma_5 \epsilon_{ijk} \gamma_j \nabla_k\f$  
   */
  class MesB1xNablaT1A0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaT1A0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (B1xNabla_T1) wavefunction onto a a1
  /*!
   * WaveFunction is  b1 x nabla_T1 onto a a1
   * The interpolator is   
   * \f$\gamma_4 \gamma_5 \epsilon_{ijk} \gamma_j \nabla_k\f$  
   */
  class MesB1xNablaT1A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaT1A1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (B1xNabla_T1) wavefunction onto a pi1
  /*!
   * WaveFunction is  b1 x nabla_T1 onto a pi1
   * The interpolator is   
   * \f$\gamma_4 \gamma_5 \epsilon_{ijk} \gamma_j \nabla_k\f$  
   */
  class MesB1xNablaT1Pi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaT1Pi1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


 //----------------------------------------------------------------------------------
  //! Construct (B1xNabla_T1) wavefunction onto a a_2
  /*!
   * WaveFunction is  b1 x nabla_T1 onto a a_2
   * The interpolator is   
   * \f$\gamma_4 \gamma_5 \epsilon_{ijk} \gamma_j \nabla_k\f$  
   */
  class MesB1xNablaT1A2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaT1A2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (B1xNabla_T1) wavefunction onto a \pi_2
  /*!
   * WaveFunction is  b1 x nabla_T1 onto a \pi_2
   * The interpolator is   
   * \f$\gamma_4 \gamma_5 \epsilon_{ijk} \gamma_j \nabla_k\f$  
   */
  class MesB1xNablaT1Pi2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaT1Pi2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };





 //----------------------------------------------------------------------------------
  //! Construct (B1xNabla_T2) wavefunction onto a a1
  /*!
   * WaveFunction is  b1 x nabla_T2 onto a a1
   * The interpolator is   
   * \f$\gamma_4 \gamma_5 s_{ijk} \gamma_j \nabla_k\f$  
   */
  class MesB1xNablaT2A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaT2A1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (B1xNabla_T2) wavefunction onto a pi1
  /*!
   * WaveFunction is  b1 x nabla_T2 onto a pi1
   * The interpolator is   
   * \f$\gamma_4 \gamma_5 s_{ijk} \gamma_j \nabla_k\f$  
   */
  class MesB1xNablaT2Pi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaT2Pi1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (B1xNabla_T2) wavefunction onto a a2
  /*!
   * WaveFunction is  b1 x nabla_T2 onto a a2
   * The interpolator is   
   * \f$\gamma_4 \gamma_5 s_{ijk} \gamma_j \nabla_k\f$  
   */
  class MesB1xNablaT2A2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaT2A2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (B1xNabla_T2) wavefunction onto a pi2
  /*!
   * WaveFunction is  b1 x nabla_T2 onto a pi2
   * The interpolator is   
   * \f$\gamma_4 \gamma_5 s_{ijk} \gamma_j \nabla_k\f$  
   */
  class MesB1xNablaT2Pi2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaT2Pi2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


 //----------------------------------------------------------------------------------
  //! Construct (B1xNabla_E) wavefunction onto a a1
  /*!
   * WaveFunction is  b1 x nabla_E onto a a1
   * The interpolator is   
   * \f$\gamma_4 \gamma_5 S_{ijk} \gamma_j \nabla_k\f$  
   */
  class MesB1xNablaEA1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaEA1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (B1xNabla_E) wavefunction onto a pi1
  /*!
   * WaveFunction is  b1 x nabla_E onto a pi1
   * The interpolator is   
   * \f$\gamma_4 \gamma_5 S_{ijk} \gamma_j \nabla_k\f$  
   */
  class MesB1xNablaEPi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaEPi1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (B1xNabla_E) wavefunction onto a a2
  /*!
   * WaveFunction is  b1 x nabla_E onto a a2
   * The interpolator is   
   * \f$\gamma_4 \gamma_5 S_{ijk} \gamma_j \nabla_k\f$  
   */
  class MesB1xNablaEA2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaEA2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (B1xNabla_E) wavefunction onto a pi2
  /*!
   * WaveFunction is  b1 x nabla_E onto a pi2
   * The interpolator is   
   * \f$\gamma_4 \gamma_5 S_{ijk} \gamma_j \nabla_k\f$  
   */
  class MesB1xNablaEPi2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xNablaEPi2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (A0_2xD_T2) wavefunction onto b_0
  /*!
   * WaveFunction is  a0_2 x D_T2 onto b_0 
   * The interpolator is   
   * \f$\gamma_4 D_i\f$  
   */
  class MesA02xDT2B0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA02xDT2B0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (A0_2xD_T2) wavefunction onto b_1
  /*!
   * WaveFunction is  a0_2 x D_T2 onto b_1
   * The interpolator is   
   * \f$\gamma_4 D_i\f$  
   */
  class MesA02xDT2B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA02xDT2B1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (A0_2xD_T2) wavefunction onto \rho_1
  /*!
   * WaveFunction is  a0_2 x D_T2 onto \rho_1
   * The interpolator is   
   * \f$\gamma_4 D_i\f$  
   */
  class MesA02xDT2Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA02xDT2Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------
  //! Construct (A0_2xD_T2) wavefunction onto b_2
  /*!
   * WaveFunction is  a0_2 x D_T2 onto b_2
   * The interpolator is   
   * \f$\gamma_4 D_i\f$  
   */
  class MesA02xDT2B2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA02xDT2B2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


 //----------------------------------------------------------------------------------
  //! Construct (A0_2xD_T2) wavefunction onto \rho_2
  /*!
   * WaveFunction is  a0_2 x D_T2 onto \rho_2
   * The interpolator is   
   * \f$\gamma_4 D_i\f$  
   */
  class MesA02xDT2Rho2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA02xDT2Rho2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


 //----------------------------------------------------------------------------------
  //! Construct (A0_2xD_T2) wavefunction onto \rho_3
  /*!
   * WaveFunction is  a0_2 x D_T2 onto \rho_3
   * The interpolator is   
   * \f$\gamma_4 D_i\f$  
   */
  class MesA02xDT2Rho3WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA02xDT2Rho3WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho3";}

    //! The number of polarizations
    int numPolar() const {return 7;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };



  //-----------------------------------------------------------------------
  //! Construct (A1xD_A2) wavefunction onto \pi_2
  /*!
   * WaveFunction is  a1 x D_A2 onto \pi_2
   * The interpolator is   
   * \f$\gamma_5\gamma_i D_i\f$  
   */
  class MesA1xDA2Pi2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDA2Pi2WaveFunc(const WaveFuncParams& p) : params(p) {}

   //! The name of this particle
    std::string particleName() const {return "pi2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //-----------------------------------------------------------------------
  //! Construct (A1xD_A2) wavefunction onto a_3
  /*!
   * WaveFunction is  a1 x D_A2 onto a_3
   * The interpolator is   
   * \f$\gamma_5\gamma_i D_i\f$  
   */
  class MesA1xDA2A3WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDA2A3WaveFunc(const WaveFuncParams& p) : params(p) {}

   //! The name of this particle
    std::string particleName() const {return "a3";}

    //! The number of polarizations
    int numPolar() const {return 7;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //--------------------------------------------------------------------------------
  //! Construct (A1xD_T1) wavefunction onto a \pi_0
  /*!
   * WaveFunction is  a1 x D_T1 onto a \pi_0
   * The interpolator is   
   * \f$\gamma_5 s_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDT1Pi0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDT1Pi0WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //--------------------------------------------------------------------------------
  //! Construct (A1xD_T1) wavefunction onto a \pi_1
  /*!
   * WaveFunction is  a1 x D_T1 onto a \pi_1
   * The interpolator is   
   * \f$\gamma_5 s_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDT1Pi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDT1Pi1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //--------------------------------------------------------------------------------
  //! Construct (A1xD_T1) wavefunction onto a a_1
  /*!
   * WaveFunction is  a1 x D_T1 onto a a_1
   * The interpolator is   
   * \f$\gamma_5 s_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDT1A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDT1A1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //--------------------------------------------------------------------------------
  //! Construct (A1xD_T1) wavefunction onto a \pi_2
  /*!
   * WaveFunction is  a1 x D_T1 onto a \pi_2
   * The interpolator is   
   * \f$\gamma_5 s_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDT1Pi2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDT1Pi2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //--------------------------------------------------------------------------------
  //! Construct (A1xD_T1) wavefunction onto a a_2
  /*!
   * WaveFunction is  a1 x D_T1 onto a a_2
   * The interpolator is   
   * \f$\gamma_5 s_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDT1A2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDT1A2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //--------------------------------------------------------------------------------
  //! Construct (A1xD_T1) wavefunction onto a a_3
  /*!
   * WaveFunction is  a1 x D_T1 onto a a_3
   * The interpolator is   
   * \f$\gamma_5 s_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDT1A3WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDT1A3WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a3";}

    //! The number of polarizations
    int numPolar() const {return 7;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //----------------------------------------------------------------------------------------
  //! Construct (A1xD_T2) wavefunction onto a \pi0
  /*!
   * WaveFunction is  a1 x D_T2 onto a \pi0
   * The interpolator is   
   * \f$\gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDT2Pi0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDT2Pi0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


 //----------------------------------------------------------------------------------------
  //! Construct (A1xD_T2) wavefunction onto a \pi1
  /*!
   * WaveFunction is  a1 x D_T2 onto a \pi1
   * The interpolator is   
   * \f$\gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDT2Pi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDT2Pi1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------------
  //! Construct (A1xD_T2) wavefunction onto a a1
  /*!
   * WaveFunction is  a1 x D_T2 onto a a1
   * The interpolator is   
   * \f$\gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDT2A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDT2A1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //----------------------------------------------------------------------------------------
  //! Construct (A1xD_T2) wavefunction onto a \pi_2
  /*!
   * WaveFunction is  a1 x D_T2 onto a \pi_2
   * The interpolator is 
   * \f$\gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDT2Pi2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDT2Pi2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


 //----------------------------------------------------------------------------------------
  //! Construct (A1xD_T2) wavefunction onto a a_2
  /*!
   * WaveFunction is  a1 x D_T2 onto a a_2
   * The interpolator is 
   * \f$\gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDT2A2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDT2A2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//----------------------------------------------------------------------------------------
  //! Construct (A1xD_T2) wavefunction onto a a_3
  /*!
   * WaveFunction is  a1 x D_T2 onto a a_3
   * The interpolator is 
   * \f$\gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDT2A3WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDT2A3WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a3";}

    //! The number of polarizations
    int numPolar() const {return 7;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


//----------------------------------------------------------------------------------------
  //! Construct (A1xD_E) wavefunction onto a \pi_1
  /*!
   * WaveFunction is  a1 x D_E onto a \pi_1
   * The interpolator is 
   * \f$\gamma_5 S_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDEPi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDEPi1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


//----------------------------------------------------------------------------------------
  //! Construct (A1xD_E) wavefunction onto a a_1
  /*!
   * WaveFunction is  a1 x D_E onto a a_1
   * The interpolator is 
   * \f$\gamma_5 S_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDEA1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDEA1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


//----------------------------------------------------------------------------------------
  //! Construct (A1xD_E) wavefunction onto a \pi_2
  /*!
   * WaveFunction is  a1 x D_E onto a \pi_2
   * The interpolator is 
   * \f$\gamma_5 S_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDEPi2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDEPi2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "pi2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//----------------------------------------------------------------------------------------
  //! Construct (A1xD_E) wavefunction onto a a_2
  /*!
   * WaveFunction is  a1 x D_E onto a a_2
   * The interpolator is 
   * \f$\gamma_5 S_{ijk}\gamma_j D_k\f$  
   */
  class MesA1xDEA2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xDEA2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "a2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };



#if 0

  //! Construct (B1xD_A2) wavefunction
  /*!
   * WaveFunction is  b1 x D_A2
   * The interpolator is   
   * \f$\gamma_4\gamma_5 \gamma_i D_i\f$  
   */
  class MesB1xDA2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xDA2WaveFunc(const WaveFuncParams& p) : params(p) {}

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //! Construct (B1xD_E) wavefunction
  /*!
   * WaveFunction is  b1 x D_E
   * The interpolator is   
   * \f$\gamma_4\gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
   */
  class MesB1xDEWaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xDEWaveFunc(const WaveFuncParams& p) : params(p) {}

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (B1xD_T1) wavefunction
  /*!
   * WaveFunction is  b1 x D_T1
   * The interpolator is   
   * \f$\gamma_4\gamma_5 s_{ijk}\gamma_j D_k\f$  
   */
  class MesB1xDT1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xDT1WaveFunc(const WaveFuncParams& p) : params(p) {}

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //! Construct (B1xD_T2) wavefunction
  /*!
   * WaveFunction is  b1 x D_T2
   * The interpolator is   
   * \f$\gamma_4\gamma_5 \epsilon_{ijk}\gamma_j D_k\f$  
   */
  class MesB1xDT2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesB1xDT2WaveFunc(const WaveFuncParams& p) : params(p) {}

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

#endif

  //-------------------------------------------------------------------------------------
  //! Construct (RhoxD_A2) wavefunction onto a b_2
  /*!
   * WaveFunction is  rho x D_A2 onto a b_2
   * The interpolator is   
   * \f$\gamma_i D_i\f$  
   */
  class MesRhoxDA2B2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDA2B2WaveFunc(const WaveFuncParams& p) : params(p) {}

   //! The name of this particle
    std::string particleName() const {return "b2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //-------------------------------------------------------------------------------------
  //! Construct (RhoxD_A2) wavefunction onto a \rho_3
  /*!
   * WaveFunction is  rho x D_A2 onto a \rho_3
   * The interpolator is   
   * \f$\gamma_i D_i\f$  
   */
  class MesRhoxDA2Rho3WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDA2Rho3WaveFunc(const WaveFuncParams& p) : params(p) {}

   //! The name of this particle
    std::string particleName() const {return "rho3";}

    //! The number of polarizations
    int numPolar() const {return 7;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };



  //---------------------------------------------------------------------------------
  //! Construct (RhoxD_T1) wavefunction onto a b_0
  /*!
   * WaveFunction is  rho x D_T1 onto a b_0 
   * The interpolator is   
   * \f$s_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDT1B0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDT1B0WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //---------------------------------------------------------------------------------
  //! Construct (RhoxD_T1) wavefunction onto a b_1
  /*!
   * WaveFunction is  rho x D_T1 onto a b_1
   * The interpolator is   
   * \f$s_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDT1B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDT1B1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //---------------------------------------------------------------------------------
  //! Construct (RhoxD_T1) wavefunction onto a \rho_1
  /*!
   * WaveFunction is  rho x D_T1 onto a \rho_1
   * The interpolator is   
   * \f$s_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDT1Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDT1Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //---------------------------------------------------------------------------------
  //! Construct (RhoxD_T1) wavefunction onto a b_2
  /*!
   * WaveFunction is  rho x D_T1 onto a b_2
   * The interpolator is   
   * \f$s_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDT1B2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDT1B2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "b2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //---------------------------------------------------------------------------------
  //! Construct (RhoxD_T1) wavefunction onto a \rho_2
  /*!
   * WaveFunction is  rho x D_T1 onto a \rho_2
   * The interpolator is   
   * \f$s_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDT1Rho2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDT1Rho2WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //---------------------------------------------------------------------------------
  //! Construct (RhoxD_T1) wavefunction onto a \rho_3
  /*!
   * WaveFunction is  rho x D_T1 onto a \rho_3
   * The interpolator is   
   * \f$s_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDT1Rho3WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDT1Rho3WaveFunc(const WaveFuncParams& p) : params(p) {}

    //! The name of this particle
    std::string particleName() const {return "rho3";}

    //! The number of polarizations
    int numPolar() const {return 7;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;

  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //------------------------------------------------------------------------------------------
  //! Construct (RhoxD_T2) wavefunction onto a b_0
  /*!
   * WaveFunction is  rho x D_T2 onto a b_0
   * The interpolator is   
   * \f$\epsilon_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDT2B0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDT2B0WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


  //------------------------------------------------------------------------------------------
  //! Construct (RhoxD_T2) wavefunction onto a b_1
  /*!
   * WaveFunction is  rho x D_T2 onto a b_1
   * The interpolator is   
   * \f$\epsilon_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDT2B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDT2B1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //------------------------------------------------------------------------------------------
  //! Construct (RhoxD_T2) wavefunction onto a \rho_1
  /*!
   * WaveFunction is  rho x D_T2 onto a \rho_1
   * The interpolator is   
   * \f$\epsilon_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDT2Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDT2Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


 //------------------------------------------------------------------------------------------
  //! Construct (RhoxD_T2) wavefunction onto a b_2
  /*!
   * WaveFunction is  rho x D_T2 onto a b_2
   * The interpolator is
   * \f$\epsilon_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDT2B2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDT2B2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //------------------------------------------------------------------------------------------
  //! Construct (RhoxD_T2) wavefunction onto a \rho_2
  /*!
   * WaveFunction is  rho x D_T2 onto a \rho_2
   * The interpolator is
   * \f$\epsilon_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDT2Rho2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDT2Rho2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //------------------------------------------------------------------------------------------
  //! Construct (RhoxD_T2) wavefunction onto a \rho_3
  /*!
   * WaveFunction is  rho x D_T2 onto a \rho_3
   * The interpolator is
   * \f$\epsilon_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDT2Rho3WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDT2Rho3WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho3";}

    //! The number of polarizations
    int numPolar() const {return 7;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };




 //------------------------------------------------------------------------------------------
  //! Construct (RhoxD_E) wavefunction onto a b_1
  /*!
   * WaveFunction is  rho x D_E onto a b_1
   * The interpolator is
   * \f$S_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDEB1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDEB1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


 //------------------------------------------------------------------------------------------
  //! Construct (RhoxD_E) wavefunction onto a rho_1
  /*!
   * WaveFunction is  rho x D_E onto a rho_1
   * The interpolator is
   * \f$S_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDERho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDERho1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //------------------------------------------------------------------------------------------
  //! Construct (RhoxD_E) wavefunction onto a b_2
  /*!
   * WaveFunction is  rho x D_E onto a b_2
   * The interpolator is
   * \f$S_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDEB2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDEB2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


 //------------------------------------------------------------------------------------------
  //! Construct (RhoxD_E) wavefunction onto a \rho_2
  /*!
   * WaveFunction is  rho x D_E onto a \rho_2
   * The interpolator is
   * \f$S_{ijk}\gamma_j D_k\f$  
   */
  class MesRhoxDERho2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxDERho2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 2;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //------------------------------------------------------------------------------------------
  //! Construct (Pion_2xD_T2) wavefunction onto a \pi_0
  /*!
   * WaveFunction is  pion2 x D_T2 onto a \pi_0
   * The interpolator is
   * \f$\gamma_4 \gamma_5 D_k\f$  
   */
  class MesPion2xDT2Pi0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesPion2xDT2Pi0WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //------------------------------------------------------------------------------------------
  //! Construct (Pion_2xD_T2) wavefunction onto a \pi_1
  /*!
   * WaveFunction is  pion2 x D_T2 onto a \pi_1
   * The interpolator is
   * \f$\gamma_4 \gamma_5 D_k\f$  
   */
  class MesPion2xDT2Pi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesPion2xDT2Pi1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //------------------------------------------------------------------------------------------
  //! Construct (Pion_2xD_T2) wavefunction onto a a_1
  /*!
   * WaveFunction is  pion2 x D_T2 onto a a_1
   * The interpolator is
   * \f$\gamma_4 \gamma_5 D_k\f$  
   */
  class MesPion2xDT2A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesPion2xDT2A1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //------------------------------------------------------------------------------------------
  //! Construct (Pion_2xD_T2) wavefunction onto a \pi_2
  /*!
   * WaveFunction is  pion2 x D_T2 onto a \pi_2
   * The interpolator is
   * \f$\gamma_4 \gamma_5 D_k\f$  
   */
  class MesPion2xDT2Pi2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesPion2xDT2Pi2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //------------------------------------------------------------------------------------------
  //! Construct (Pion_2xD_T2) wavefunction onto a a_2
  /*!
   * WaveFunction is  pion2 x D_T2 onto a a_2
   * The interpolator is
   * \f$\gamma_4 \gamma_5 D_k\f$  
   */
  class MesPion2xDT2A2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesPion2xDT2A2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //------------------------------------------------------------------------------------------
  //! Construct (Pion_2xD_T2) wavefunction onto a a_3
  /*!
   * WaveFunction is  pion2 x D_T2 onto a a_3
   * The interpolator is
   * \f$\gamma_4 \gamma_5 D_k\f$  
   */
  class MesPion2xDT2A3WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesPion2xDT2A3WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a3";}

    //! The number of polarizations
    int numPolar() const {return 7;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //---------------------------------------------------------------------------------
  //! Construct (PionxD_T2) wavefunction onto a \pi_0
  /*!
   * WaveFunction is  pion x D_T2 onto a \pi_0
   * The interpolator is
   * \f$\gamma_5 D_k\f$  
   */
  class MesPionxDT2Pi0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesPionxDT2Pi0WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //---------------------------------------------------------------------------------
  //! Construct (PionxD_T2) wavefunction onto a \pi_2
  /*!
   * WaveFunction is  pion x D_T2 onto a \pi_2
   * The interpolator is
   * \f$\gamma_5 D_k\f$  
   */
  class MesPionxDT2Pi2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesPionxDT2Pi2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //---------------------------------------------------------------------------------
  //! Construct (PionxD_T2) wavefunction onto a a_1
  /*!
   * WaveFunction is  pion x D_T2 onto a a_1
   * The interpolator is
   * \f$\gamma_5 D_k\f$  
   */
  class MesPionxDT2A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesPionxDT2A1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //---------------------------------------------------------------------------------
  //! Construct (A0xD_T2) wavefunction onto a a_0
  /*!
   * WaveFunction is  a0 x D_T2 onto a a_0
   * The interpolator is
   * \f$D_k\f$  
   */
  class MesA0xDT2A0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA0xDT2A0WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //---------------------------------------------------------------------------------
  //! Construct (A0xD_T2) wavefunction onto a a_2
  /*!
   * WaveFunction is  a0 x D_T2 onto a a_2
   * The interpolator is
   * \f$D_k\f$  
   */
  class MesA0xDT2A2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA0xDT2A2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //---------------------------------------------------------------------------------
  //! Construct (A0xD_T2) wavefunction onto a pi_1
  /*!
   * WaveFunction is  a0 x D_T2 onto a pi_1
   * The interpolator is
   * \f$D_k\f$  
   */
  class MesA0xDT2Pi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA0xDT2Pi1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

  //---------------------------------------------------------------------------------
  //! Construct (A0xB_T1) wavefunction onto a b_1
  /*!
   * WaveFunction is  a0 x B_T1 onto a b_1
   * The interpolator is
   * \f$B_k\f$  
   */
  class MesA0xBT1B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA0xBT1B1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //---------------------------------------------------------------------------------
  //! Construct (A0xB_T1) wavefunction onto a \rho_1
  /*!
   * WaveFunction is  a0 x B_T1 onto a \rho_1
   * The interpolator is
   * \f$B_k\f$  
   */
  class MesA0xBT1Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA0xBT1Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //---------------------------------------------------------------------------------
  //! Construct (PionxB_T1) wavefunction onto a \rho_1
  /*!
   * WaveFunction is  a0 x B_T1 onto a \rho_1
   * The interpolator is
   * \f$\gamma_5 B_k\f$  
   */
  class MesPionxBT1Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesPionxBT1Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //---------------------------------------------------------------------------------
  //! Construct (PionxB_T1) wavefunction onto a b_1
  /*!
   * WaveFunction is  a0 x B_T1 onto a b_1
   * The interpolator is
   * \f$\gamma_5 B_k\f$  
   */
  class MesPionxBT1B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesPionxBT1B1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //---------------------------------------------------------------------------------
  //! Construct (A0_2xB_T1) wavefunction onto a \pi_0
  /*!
   * WaveFunction is  a02 x B_T1 onto a \pi_0
   * The interpolator is
   * \f$\gamma_5 \gamma_4 B_k\f$  
   */
  class MesA02xBT1Pi0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA02xBT1Pi0WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //---------------------------------------------------------------------------------
  //! Construct (A0_2xB_T1) wavefunction onto a a_1
  /*!
   * WaveFunction is  a02 x B_T1 onto a a_1
   * The interpolator is
   * \f$\gamma_5 \gamma_4 B_k\f$  
   */
  class MesA02xBT1A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA02xBT1A1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


 //---------------------------------------------------------------------------------
  //! Construct (A0_2xB_T1) wavefunction onto a \pi_1
  /*!
   * WaveFunction is  a02 x B_T1 onto a \pi_1
   * The interpolator is
   * \f$\gamma_5 \gamma_4 B_k\f$  
   */
  class MesA02xBT1Pi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA02xBT1Pi1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


 //---------------------------------------------------------------------------------
  //! Construct (A0_2xB_T1) wavefunction onto a a_2
  /*!
   * WaveFunction is  a02 x B_T1 onto a a_2
   * The interpolator is
   * \f$\gamma_5 \gamma_4 B_k\f$  
   */
  class MesA02xBT1A2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA02xBT1A2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //---------------------------------------------------------------------------------
  //! Construct (A0_2xB_T1) wavefunction onto a \pi_2
  /*!
   * WaveFunction is  a02 x B_T1 onto a \pi_2
   * The interpolator is
   * \f$\gamma_5 \gamma_4 B_k\f$  
   */
  class MesA02xBT1Pi2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA02xBT1Pi2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


 //---------------------------------------------------------------------------------
  //! Construct (RhoxB_A1) wavefunction onto a \pi_0
  /*!
   * WaveFunction is  rho x B_A1 onto a \pi_0
   * The interpolator is
   * \f$\gamma_k B_k\f$  
   */
  class MesRhoxBA1Pi0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBA1Pi0WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //---------------------------------------------------------------------------------
  //! Construct (RhoxB_A1) wavefunction onto a a_1
  /*!
   * WaveFunction is  rho x B_A1 onto a a_1
   * The interpolator is
   * \f$\gamma_k B_k\f$  
   */
  class MesRhoxBA1A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBA1A1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //---------------------------------------------------------------------------------
  //! Construct (RhoxB_A1) wavefunction onto a \pi_2
  /*!
   * WaveFunction is  rho x B_A1 onto a \pi_2
   * The interpolator is
   * \f$\gamma_k B_k\f$  
   */
  class MesRhoxBA1Pi2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBA1Pi2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

 //---------------------------------------------------------------------------------
  //! Construct (RhoxB_T1) wavefunction onto a a_0
  /*!
   * WaveFunction is  rho x B_T1 onto a a_0
   * The interpolator is
   * \f$\epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesRhoxBT1A0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBT1A0WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (RhoxB_T1) wavefunction onto a a_1
  /*!
   * WaveFunction is  rho x B_T1 onto a a_1
   * The interpolator is
   * \f$\epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesRhoxBT1A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBT1A1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (RhoxB_T1) wavefunction onto a \pi_1
  /*!
   * WaveFunction is  rho x B_T1 onto a \pi_1
   * The interpolator is
   * \f$\epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesRhoxBT1Pi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBT1Pi1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (RhoxB_T1) wavefunction onto a a_2
  /*!
   * WaveFunction is  rho x B_T1 onto a a_2
   * The interpolator is
   * \f$\epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesRhoxBT1A2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBT1A2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (RhoxB_T1) wavefunction onto a \pi_2
  /*!
   * WaveFunction is  rho x B_T1 onto a \pi_2
   * The interpolator is
   * \f$\epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesRhoxBT1Pi2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBT1Pi2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (RhoxB_T2) wavefunction onto a a_1
  /*!
   * WaveFunction is  rho x B_T2 onto a a_1
   * The interpolator is
   * \f$s_{ijk} \gamma_j B_k\f$  
   */
  class MesRhoxBT2A1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBT2A1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (RhoxB_T2) wavefunction onto a \pi_1
  /*!
   * WaveFunction is  rho x B_T2 onto a \pi_1
   * The interpolator is
   * \f$s_{ijk} \gamma_j B_k\f$  
   */
  class MesRhoxBT2Pi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBT2Pi1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (RhoxB_T2) wavefunction onto a a_2
  /*!
   * WaveFunction is  rho x B_T2 onto a a_2
   * The interpolator is
   * \f$s_{ijk} \gamma_j B_k\f$  
   */
  class MesRhoxBT2A2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBT2A2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (RhoxB_T2) wavefunction onto a \pi_2
  /*!
   * WaveFunction is  rho x B_T2 onto a \pi_2
   * The interpolator is
   * \f$s_{ijk} \gamma_j B_k\f$  
   */
  class MesRhoxBT2Pi2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBT2Pi2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (RhoxB_E) wavefunction onto a a_1
  /*!
   * WaveFunction is  rho x B_E onto a a_1
   * The interpolator is
   * \f$s_{ijk} \gamma_j B_k\f$  
   */
  class MesRhoxBEA1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBEA1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (RhoxB_E) wavefunction onto a \pi_1
  /*!
   * WaveFunction is  rho x B_E onto a \pi_1
   * The interpolator is
   * \f$s_{ijk} \gamma_j B_k\f$  
   */
  class MesRhoxBEPi1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBEPi1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (RhoxB_E) wavefunction onto a a_2
  /*!
   * WaveFunction is  rho x B_E onto a a_2
   * The interpolator is
   * \f$s_{ijk} \gamma_j B_k\f$  
   */
  class MesRhoxBEA2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBEA2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "a2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (RhoxB_E) wavefunction onto a \pi_2
  /*!
   * WaveFunction is  rho x B_E onto a \pi_2
   * The interpolator is
   * \f$s_{ijk} \gamma_j B_k\f$  
   */
  class MesRhoxBEPi2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesRhoxBEPi2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "pi2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (A1xB_A1) wavefunction onto a b_0
  /*!
   * WaveFunction is  a1 x B_A1 onto a b_0
   * The interpolator is
   * \f$\gamma_5 \gamma_k B_k\f$  
   */
  class MesA1xBA1B0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xBA1B0WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (A1xB_A1) wavefunction onto a \rho_1
  /*!
   * WaveFunction is  a1 x B_A1 onto a \rho_1
   * The interpolator is
   * \f$\gamma_5 \gamma_k B_k\f$  
   */
  class MesA1xBA1Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xBA1Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (A1xB_A1) wavefunction onto a b_2
  /*!
   * WaveFunction is  a1 x B_A1 onto a b_2
   * The interpolator is
   * \f$\gamma_5 \gamma_k B_k\f$  
   */
  class MesA1xBA1B2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xBA1B2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 0;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (A1xB_T1) wavefunction onto a \rho_0
  /*!
   * WaveFunction is  a1 x B_A1 onto a \rho_0
   * The interpolator is
   * \f$\gamma_5 \epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesA1xBT1Rho0WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xBT1Rho0WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho0";}

    //! The number of polarizations
    int numPolar() const {return 1;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (A1xB_T1) wavefunction onto a \rho_1
  /*!
   * WaveFunction is  a1 x B_A1 onto a \rho_1
   * The interpolator is
   * \f$\gamma_5 \epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesA1xBT1Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xBT1Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (A1xB_T1) wavefunction onto a b_1
  /*!
   * WaveFunction is  a1 x B_A1 onto a b_1
   * The interpolator is
   * \f$\gamma_5 \epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesA1xBT1B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xBT1B1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (A1xB_T1) wavefunction onto a \rho_2
  /*!
   * WaveFunction is  a1 x B_A1 onto a \rho_2
   * The interpolator is
   * \f$\gamma_5 \epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesA1xBT1Rho2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xBT1Rho2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (A1xB_T1) wavefunction onto a b_2
  /*!
   * WaveFunction is  a1 x B_A1 onto a b_2
   * The interpolator is
   * \f$\gamma_5 \epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesA1xBT1B2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xBT1B2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (A1xB_T2) wavefunction onto a \rho_1
  /*!
   * WaveFunction is  a1 x B_A1 onto a \rho_1
   * The interpolator is
   * \f$\gamma_5 \epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesA1xBT2Rho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xBT2Rho1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (A1xB_T2) wavefunction onto a b_1
  /*!
   * WaveFunction is  a1 x B_A1 onto a b_1
   * The interpolator is
   * \f$\gamma_5 \epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesA1xBT2B1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xBT2B1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (A1xB_T2) wavefunction onto a \rho_2
  /*!
   * WaveFunction is  a1 x B_A1 onto a \rho_2
   * The interpolator is
   * \f$\gamma_5 \epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesA1xBT2Rho2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xBT2Rho2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (A1xB_T2) wavefunction onto a b_2
  /*!
   * WaveFunction is  a1 x B_A1 onto a b_2
   * The interpolator is
   * \f$\gamma_5 \epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesA1xBT2B2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xBT2B2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (A1xB_E) wavefunction onto a \rho_1
  /*!
   * WaveFunction is  a1 x B_A1 onto a \rho_1
   * The interpolator is
   * \f$\gamma_5 \epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesA1xBERho1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xBERho1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (A1xB_E) wavefunction onto a b_1
  /*!
   * WaveFunction is  a1 x B_A1 onto a b_1
   * The interpolator is
   * \f$\gamma_5 \epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesA1xBEB1WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xBEB1WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b1";}

    //! The number of polarizations
    int numPolar() const {return 3;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (A1xB_E) wavefunction onto a \rho_2
  /*!
   * WaveFunction is  a1 x B_A1 onto a \rho_2
   * The interpolator is
   * \f$\gamma_5 \epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesA1xBERho2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xBERho2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "rho2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };

//---------------------------------------------------------------------------------
  //! Construct (A1xB_E) wavefunction onto a b_2
  /*!
   * WaveFunction is  a1 x B_A1 onto a b_2
   * The interpolator is
   * \f$\gamma_5 \epsilon_{ijk} \gamma_j B_k\f$  
   */
  class MesA1xBEB2WaveFunc : public WaveFunction<Complex>
  {
  public:
    //! Full constructor
    MesA1xBEB2WaveFunc(const WaveFuncParams& p) : params(p) {}

  //! The name of this particle
    std::string particleName() const {return "b2";}

    //! The number of polarizations
    int numPolar() const {return 5;}

    //! The number of directions
    int numDir() const {return 3;}

    //! The overlap factor
    Array<Complex> operator()(const Real& mass, const ArrayInt& p, int dir, int r) const;
 
  private:
    WaveFuncParams  params;   /*!< wavefunction params */
  };


} // namespace FF

#endif
