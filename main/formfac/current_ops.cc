// $Id: current_ops.cc,v 2.0 2008/12/05 04:43:47 edwards Exp $
/*! \file
 * \brief Current operators
 */

#include "current_ops.h"
#include "current_ops_factory.h"
#include "parton/parton_distribution_moments.h"
#include "formfac/formfac_manage_factory.h"
#include "formfac/formfac_manage_aggregate.h"

namespace FF
{
  //----------------------------------------------------------------------------
  //! Construct a prototype key for 3-pts
  ThreePtArg constructPrototype(const LLSqRow_t& in)
  {
    ThreePtArg out;

    out.dt          = in.dt;
    out.quark       = in.quark;

    out.src_name    = in.src_name;
    out.src_smear   = in.src_smear;
    out.src_lorentz = in.src_lorentz;
    out.src_spin    = in.src_spin;

    out.snk_name    = in.snk_name;
    out.snk_smear   = in.snk_smear;
    out.snk_lorentz = in.snk_lorentz;
    out.snk_spin    = in.snk_spin;

    out.mass        = in.mass;
    out.ensemble    = in.ensemble;

    return out;
  }


  //---------------------------------------------------------------------------
  //! Convenience - make a reduced threept arg
  ThreePtArgReduced threePtArg(const PiPf& pi_pf, int g)
  {
    ThreePtArgReduced  arg;
    arg.pi_pf  = pi_pf;
    arg.gamma  = g; 

    return arg;
  }

  //! Convenience - make a reduced threept arg with link insertions
  ThreePtArgReduced threePtArg(const PiPf& pi_pf, int g, int link_dir)
  {
    ThreePtArgReduced  arg;
    arg.pi_pf  = pi_pf;
    arg.gamma  = g;

    arg.links.resize(1);
    arg.links[0] = link_dir;

    return arg;
  }

  //! Convenience - make a reduced threept arg with link insertions
  ThreePtArgReduced threePtArg(const PiPf& pi_pf, int g, const Array<int>& links)
  {
    ThreePtArgReduced  arg;
    arg.pi_pf  = pi_pf;
    arg.gamma  = g;
    arg.links  = links;

    return arg;
  }


  //---------------------------------------------------------------------------
  //! Local vector current
  EnsemReal local_vector(Manage3PtFuncReduced& threept, const PiPf& pi_pf, int mu, int t)
  {
    // NOTE: Need to worry about signs and imag and real here 
    EnsemVectorComplex three = threept[threePtArg(pi_pf, 1 << mu)];

    // Extract appropriate part of 3pt.
    EnsemReal thr;
    if (mu == Nd-1)
      thr =  real(peekObs(three,t));
    else
      thr =  imag(peekObs(three,t));

    return thr;
  }


#if 0
// THIS STUFF NEEDS WORK
  EnsemReal conserved_vector(Manage3PtFuncReduced& threept, const PiPf& pi_pf, int mu, int t)
  {
    //constructs the symmetrised conserved current on an anisotropic lattice
    EnsemVectorComplexF three;

    if( mu == Nd-1){ //TEMPORAL
      EnsemVectorComplexF three_a = project( threePtArg(pi_pf[p], 0,  nu, rho, 3), p);
      EnsemVectorComplexF three_b = project( threePtArg(pi_pf[p], 8,  nu, rho, 3), p);
      EnsemVectorComplexF three_c = project( threePtArg(pi_pf[p], 0,  nu, rho, 7), p);
      EnsemVectorComplexF three_d = project( threePtArg(pi_pf[p], 8,  nu, rho, 7), p);
    
      three = Real(0.25) * ( shift( three_a + three_b , 1 ) + (three_a + three_b) 
			     - ( three_c - three_d ) - shift(three_c - three_d, -1)  );
    }
    else{ //SPATIAL
      EnsemVectorComplexF three_a = project( threePtArg(pi_pf[p], 0, nu, rho, mu ), p) ;
      EnsemVectorComplexF three_b = project( threePtArg(pi_pf[p], 1 << mu, nu, rho, mu ), p);
      EnsemVectorComplexF three_c = project( threePtArg(pi_pf[p], 0, nu, rho, mu + 4 ), p);
      EnsemVectorComplexF three_d = project( threePtArg(pi_pf[p], 1 << mu, nu, rho, mu + 4 ), p);
      
      const double xi  = param.lattice.aniso.xi;
      const double csq = param.lattice.aniso.c_sq;
      const double xi_0 = param.lattice.aniso.xi_0;
      const double nu_ferm = param.lattice.aniso.nu;

      //phase is so that q = p_i - p_f;
      ArrayDouble pp_i = sqrt(csq) * contMom(pi_pf[p].p_i, param.lattice.latt_size[0]); 
      ArrayDouble pp_f = sqrt(csq) * contMom(pi_pf[p].p_f, param.lattice.latt_size[0]);
      ArrayDouble qq = pp_i - pp_f;
      Complex phase_plus = cmplx(  cos( Real(qq[mu]) ), sin( Real(qq[mu]) )  );
      Complex phase_minus = cmplx(   cos( Real(qq[mu]) ), Real(-1.0)*sin( Real(qq[mu]) )  );
      Complex one = cmplx( Real(1.0) , Real(0.0) );

      three = Real(0.25) * ( Real(nu_ferm) * Real(xi) / Real(xi_0) ) * ( (one + phase_minus )*(three_a + three_b) -  (one + phase_plus ) * ( three_c - three_d ) );    
    }
    return three;
  }

  std::vector<EnsemVectorComplexF> improvement_terms(int mu, int nu, int rho, int p)
  {
    // implements JJD's understanding of the Kronfeld, Harada... anisotropic improvement with r=1
    // NB this is not space-time symmetric obviously
    // and it doesn't claim to be any kind of conserved current

    // in this implementation i've liberally used (tree-level) equations of motion to simplify the improvement
    // it only uses the tensor terms, not requiring any links

    // V^{imp}_4 = ( \bar{\psi} \gamma_4 \psi - d1 \sum_j as\partial_j \bar{\psi} \sigma_{j4} \psi )
    // V^{imp}_i = ( \bar{\psi} \gamma_4 \psi - d1 (xi0/nu) at\partial_4 \bar{\psi} \sigma_{i4} \psi )

    // there is an extra (1 + m0at) factor that is supposed to correctly (re)normalise the currents, i.e. to ensure Z_V=1 at tree level
    // in practice we'll always fix this up non-perturbatively using the pion ff at zero momentum

  
    EnsemVectorComplexF three;
    std::vector< EnsemVectorComplexF > improvements;
  
    ArrayInt sigma_rho_4(3); sigma_rho_4[0] = 9; sigma_rho_4[1] = 10; sigma_rho_4[2] = 12;

    const double csq = param.lattice.aniso.c_sq;
    //phase is so that q = p_i - p_f   - disagrees with most of JJD's notes and comments currently written into bb inline code
    ArrayDouble pp_i = sqrt(csq) * contMom(pi_pf[p].p_i, param.lattice.latt_size[0]); 
    ArrayDouble pp_f = sqrt(csq) * contMom(pi_pf[p].p_f, param.lattice.latt_size[0]);
    ArrayDouble qq = pp_i - pp_f;

    Complex one = cmplx( Real(1.0) , Real(0.0) );

    //tree level improvement coefficients
    const double m0_at = param.lattice.aniso.m0_at;
    const double xi  = param.lattice.aniso.xi;
    const double xi_0  = param.lattice.aniso.xi_0;
    const double nu_ferm  = param.lattice.aniso.nu;

    double d1 = 0.25*(1.0 - xi) + (xi*xi * m0_at / 16.0); //equation (2.33) in Kronfeld/Harada paper, n.b. only lowest order in m
    //full, all orders in m0at, result can differ if the non-pert tuned nu is used

    //  std::cout << "d1 = " << d1 << std::endl;


    if(mu == Nd-1){ //Temporal current
      three =  Real(0.0) * project( threePtArg(pi_pf[p], 8, nu, rho) , p);
      //d sigma term
      for(int j = 0; j < 3; j++){
	if( qq[j] != 0 ){
	  Complex isinq = cmplx(  Real(0.0)  , sin( Real(qq[j]) )  );
	  three += isinq * project( threePtArg( pi_pf[p], abs(sigma_rho_4[j]) , nu, rho) , p);
	}
      }
      improvements.push_back( Real(d1) * three );
    }
    else{ //Spatial current
      //d sigma term
      three = Real(0.5) * (shift(  project( threePtArg( pi_pf[p], abs(sigma_rho_4[mu]) , nu, rho) , p) , 1 ) - shift(  project( threePtArg( pi_pf[p], abs(sigma_rho_4[mu]) , nu, rho) , p), -1) );
      three *= Real( -1.0 * xi_0 / nu_ferm );
      improvements.push_back( Real(d1) * three );    
    }

    return improvements;
  }
#endif



  //---------------------------------------------------------------------------
  // Read parameters
  void read(XMLReader& xml, const std::string& path, CurrentOperatorParams& param)
  {
    CurrentOperatorParams tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const std::string& path, const CurrentOperatorParams& param)
  {
    param.writeXML(xml, path);
  }


  //! Initialize
  CurrentOperatorParams::CurrentOperatorParams()
  {
  }


  //! Read parameters
  CurrentOperatorParams::CurrentOperatorParams(XMLReader& xml, const std::string& path)
  {
    XMLReader paramtop(xml, path);

//    read(paramtop, "CurrentType", current_type);
    threept = ADATXML::readXMLGroup(paramtop, "ThreePtFuncs", "ThreePtType");
    read(paramtop, "LatticeParam", lattice);
  }

  //! Write parameters
  void CurrentOperatorParams::writeXML(XMLWriter& xml, const std::string& path) const
  {
    push(xml, path);

//      write(xml, "CurrentOperatorType",  operator_type);
    write(xml, "LatticeParam", lattice);

    pop(xml);
  }



  //---------------------------------------------------------------------------
  // Initialize a 3-pt manage object
  Manage3PtFunc* CurrentOperator::initialize3Pt(const GroupXML_t& op_xml) const
  {
    Manage3PtFunc* threept;

    // Construct 3pt manage object
    try
    {
      std::istringstream  xml_s(op_xml.xml);
      XMLReader  optop(xml_s);
	
      threept = TheManage3PtFuncFactory::Instance().createObject(op_xml.id,
								 optop,
								 op_xml.path);
    }
    catch(const std::string& e) 
    {
      std::cerr << __func__ << ": Caught Exception reading input 3pt XML: " << e << std::endl;
      exit(1);
    }

    return threept;
  }

  //---------------------------------------------------------------------------
  // Full constructor
  LocalVectorCurrentOperator::LocalVectorCurrentOperator(const CurrentOperatorParams& p) 
    : params(p), threept(initialize3Pt(params.threept)) 
  {}


  // The value at time-slice t
  EnsemReal
  LocalVectorCurrentOperator::operator()(const LLSqRow_t& row) const
  {
    if (row.insert_lorentz.size() != 1)
    {
      std::cerr << "LocalVectorCurrent: unexpected insertion array size= " 
		<< row.insert_lorentz.size() << std::endl;
      exit(1);
    }


    Handle<Manage3PtFuncReduced> threept_red(threept->createView(constructPrototype(row)));

    return local_vector(*threept_red, row.pi_pf, row.insert_lorentz[0], row.t);
  }



  //---------------------------------------------------------------------------
  // Full constructor
  ConservedVectorCurrentOperator::ConservedVectorCurrentOperator(const CurrentOperatorParams& p) 
    : params(p), threept(initialize3Pt(params.threept)) 
  {}


  // The value at time-slice t
  EnsemReal
  ConservedVectorCurrentOperator::operator()(const LLSqRow_t& row) const
  {
    if (row.insert_lorentz.size() != 1)
    {
      std::cerr << "ConservedVectorCurrent: unexpected insertion array size= " 
		<< row.insert_lorentz.size() << std::endl;
      exit(1);
    }

    std::cerr << "ConservedVecCur not implemented: also needs more params\n";
    exit(1);

    Handle<Manage3PtFuncReduced> threept_red(threept->createView(constructPrototype(row)));

    return local_vector(*threept_red, row.pi_pf, row.insert_lorentz[0], row.t);
  }



  //---------------------------------------------------------------------------
  // Full constructor
  QUpol0CurrentOperator::QUpol0CurrentOperator(const CurrentOperatorParams& p) 
    : params(p), threept(initialize3Pt(params.threept)) 
  {}


  // The value at time-slice t
  EnsemReal
  QUpol0CurrentOperator::operator()(const LLSqRow_t& row) const
  {
    if (row.insert_lorentz.size() != 0)
    {
      std::cerr << "QUpol0: unexpected insertion array size= " 
		<< row.insert_lorentz.size() << std::endl;
      exit(1);
    }

    Handle<Manage3PtFuncReduced> threept_red(threept->createView(constructPrototype(row)));

    EnsemVectorReal thr = q_upol_0(*threept_red, row.pi_pf);

    return peekObs(thr, row.t);
  }



  //! Meson version of current operators
  namespace CurrentOperatorEnv
  { 
    //! Anonymous namespace
    namespace
    {
      bool registered = false;

      //-------------------- callback functions ---------------------------------------

      //! Local current
      CurrentOperator* mesLocalVectorCurrent(XMLReader& xml_in,
					     const std::string& path)
      {
	return new LocalVectorCurrentOperator(CurrentOperatorParams(xml_in, path));
      }

      //! Conserved current
      CurrentOperator* mesConservedVectorCurrent(XMLReader& xml_in,
						 const std::string& path)
      {
	return new ConservedVectorCurrentOperator(CurrentOperatorParams(xml_in, path));
      }

      //! Upol0
      CurrentOperator* mesQUpol0Current(XMLReader& xml_in,
					const std::string& path)
      {
	return new QUpol0CurrentOperator(CurrentOperatorParams(xml_in, path));
      }

    }  // anonymous namespace


    // Register all the possible matrix elements
    bool registerAll(void) 
    {
      bool success = true;

      if (! registered)
      {
	//! Register all the factories
	success &= FormfacManage3PtFuncEnv::registerAll();

	//! Register all the new objects
	success &= TheCurrentOperatorFactory::Instance().registerObject(std::string("LocalVectorCurrent"),
									mesLocalVectorCurrent);
	success &= TheCurrentOperatorFactory::Instance().registerObject(std::string("ConservedVectorCurrent"),
									mesConservedVectorCurrent);
	success &= TheCurrentOperatorFactory::Instance().registerObject(std::string("UPol0"),
									mesLocalVectorCurrent);

	registered = true;
      }

      return success;
    }

  }  // end namespace CurrentOperatorEnv



} // namespace FF
