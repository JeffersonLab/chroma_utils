// $Id: simple_meson_2pt_w.cc,v 2.0 2008/12/05 04:44:03 edwards Exp $
/*! \file
 *  \brief Construct meson 2pt correlators.
 */

#include "simple_meson_2pt_w.h"
#include "hadron_contract_factory.h"

namespace Strippers
{

  // Read parameters
  void read(XMLReader& xml, const string& path, SimpleMeson2PtEnv::Params& param)
  {
    SimpleMeson2PtEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const SimpleMeson2PtEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Meson correlators
  /*! 
   * \ingroup hadron 
   *
   * @{
   */
  namespace SimpleMeson2PtEnv
  { 
    //! Anonymous namespace
    namespace
    {
      //-------------------- callback functions ---------------------------------------
      //! Construct pion correlator
      HadronContract* mesDiagGammaCorrs(XMLReader& xml_in,
					const std::string& path,
					int nbin)
      {
	return new DiagGammaMesonCorrs(Params(xml_in, path), nbin);   // all gammas
      }


      //! Local registration flag
      bool registered = false;

    } // end anonymous namespace


    //! Initialize
    Params::Params()
    {
    }


    //! Read parameters
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "mom2_max", mom2_max);
      read(paramtop, "avg_equiv_mom", avg_equiv_mom);
      read(paramtop, "mom_origin", mom_origin);
      read(paramtop, "first_id", first_id);
      read(paramtop, "second_id", second_id);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);

      write(xml, "mom2_max", mom2_max);
      write(xml, "avg_equiv_mom", avg_equiv_mom);
      write(xml, "mom_origin", mom_origin);

      write(xml, "first_id", first_id);
      write(xml, "second_id", second_id);

      pop(xml);
    }


    // Full constructor
    DiagGammaMesonCorrs(const Params& p, int nb) : params(p), nbin(nb)
    {
      
    }


    // Construct all the correlators
    void
    DiagGammaMesonCorrs::operator()(Handle<HadronContractResult_t> had_cont, int ibin)
    {
      QDPIO::cout << "Hadron2Pt: diagonal_gamma_mesons" << endl;

      Array<ForwardProp_t> forward_headers(2);
      forward_headers[0] = readForwardPropHeader(params.first_id);
      forward_headers[1] = readForwardPropHeader(params.second_id);
      
      Array<int> t_srce = getTSrce(forward_headers);
      int decay_dir       = getDecayDir(forward_headers);

      // Get references to the props
      const LatticePropagator& quark_prop1 = 
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.first_id);
      const LatticePropagator& quark_prop2 = 
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.second_id);

      // Parameters needed for the momentum projection
      SftMomParams_t sft_params;
      sft_params.mom2_max      = params.mom2_max;
      sft_params.origin_offset = t_srce;
      sft_params.mom_offset    = params.mom_origin;
      sft_params.avg_equiv_mom = params.avg_equiv_mom;
      sft_params.decay_dir     = decay_dir;

      std::list< Handle<Hadron2PtContract_t> > hadron;   // holds the contract lattice correlator

      for(int gamma_value=0; gamma_value < Ns*Ns; ++gamma_value)
      {
	Handle<Hadron2PtContract_t> had(new Hadron2PtContract_t);

	push(had->xml, xml_group);
	write(had->xml, id_tag, "diagonal_gamma_mesons");
	write(had->xml, "gamma_value", gamma_value);
	write(had->xml, "PropHeaders", forward_headers);
	pop(had->xml);

	had->corr = mesXCorr(quark_prop1, quark_prop2, gamma_value);

	hadron.push_back(had);  // push onto end of list
      }

      return this->project(hadron, sft_params);
    }


    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	//! Register all the factories
	success &= Strippers::TheHadronContractFactory::Instance().registerObject(string("diagonal_gamma_mesons"),
										  mesDiagGammaCorrs);

	registered = true;
      }
      return success;
    }

  }  // end namespace SimpleMeson2PtEnv

  /*! @} */   // end of group io

}  // end namespace Strippers


  
