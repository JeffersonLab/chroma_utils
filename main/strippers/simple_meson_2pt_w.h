// -*- C++ -*-
// $Id: simple_meson_2pt_w.h,v 2.0 2008/12/05 04:44:03 edwards Exp $
/*! \file
 *  \brief Read meson 2pt correlators.
 */

#ifndef __simple_meson_2pt_w_h__
#define __simple_meson_2pt_w_h__

#include "hadron_2pt.h"

namespace Strippers
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace SimpleMeson2PtEnv
  {
    bool registerAll();

    //! 2pt ensemble 
    struct Hadron2PtEnsem_t
    {
      //! Momentum projected correlator
      struct Mom_t
      {
	Array<int>          mom;    /*!< D-1 momentum of this correlator*/
	EnsemVectorComplex  corr;   /*!< Momentum projected correlator */
      };

      std::string    xml;      /*!< XML about each corr group. NOTE: THIS SHOULD DISAPPEAR */
      Array<Mom_t>   corrs;    /*!< Holds momentum projected correlators */
    };
  

    //! Simple meson 2pt parameters
    /*! @ingroup hadron */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      int           gamma_value;        /*<! Which meson is this. We know this class only produces one. */
      int           mom2_max;           /*!< (mom - mom_origin)^2 <= mom2_max */
      Array<int>    mom_origin;         /*!< Origin for the momentum */
      bool          avg_equiv_mom;      /*!< average over equivalent momenta */
    };


    //! Simple meson 2pt construction - all simple mesons
    /*! @ingroup hadron
     *
     * Create all the 16 (Ns=4) simple diagonal meson 2pt correlators
     */
    class DiagGammaMesonCorrs : public Hadron2PtCorr
    {
    public:
      //! Full constructor
      DiagGammaMesonCorrs(const Params& p, int nb);

      //! Default destructor
      ~DiagGammaMesonCorrs() {}
      
      //! Construct the correlators
      void operator()(Handle<HadronContractResult_t> had_cont, int ibin);

    private:
      //! Hide partial constructor
      DiagGammaMesonCorrs() {}

    private:
      Params   params;           /*!< The common params */
      int      nbin;             /*!< Number of jackknife bins */
      Hadron2PtEnsem_t data;     /*!< Hadron 2pt data */
    };

  }  // end namespace


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const std::string& path, SimpleMeson2PtEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const std::string& path, const SimpleMeson2PtEnv::Params& param);


}  // end namespace Strippers

#endif
