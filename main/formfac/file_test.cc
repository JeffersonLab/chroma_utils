// $Id: file_test.cc,v 2.0 2008/12/05 04:43:48 edwards Exp $

#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "wavefuncs.h"
#include "operators_factory.h"
#include <vector>
#include <iostream>

using namespace ENSEM;
using namespace FF;
using namespace std;


//! Params for linear system
struct FileTestParam
{
  FileTestParam();
  FileTestParam(XMLReader& xml_in, const std::string& path);

  LatticeParam       lattice;

  int mom2_max;
  double tol;
  string mass_string;
  string source_smear;
  string sink_smear;
};


// File test params
FileTestParam::FileTestParam(XMLReader& xml_in, const std::string& path) 
{
  try 
  {
    XMLReader paramtop(xml_in, path);

    read(paramtop, "LatticeParam", lattice);
    read(paramtop, "mom2_max", mom2_max);
    read(paramtop, "tol", tol);
    read(paramtop, "mass_string", mass_string);
    read(paramtop, "source_smear", source_smear);
    read(paramtop, "sink_smear", sink_smear);
  }
  catch(const std::string& e) 
  {
    cerr << "Caught Exception reading XML: " << e << endl;
    exit(1);
  }
}



int max(int a, int b) {return (a < b) ? b : a;}
int min(int a, int b) {return (a < b) ? a : b;}


Real mod2(const Complex& a)
{
  return real(a)*real(a) + imag(a)*imag(a);
}


string mom_name(const Array<int>& sink_mom)
{
  ostringstream s;

  if (! (sink_mom[0] == 0 && sink_mom[1] == 0 && sink_mom[2] == 0) )
    s << "_px" << sink_mom[0] << "_py" << sink_mom[1] << "_pz" << sink_mom[2];

  return s.str();
}


int main(int argc, char *argv[])
{
  cout << "Test file combinations of an operator to a gamma based operator\n\n";

  if (argc != 3)
  {
    cerr << "Usage:"
	 << argv[0] << ":  <input xml file>  <output xml file>" << endl;
    exit(1);
  }

  bool registered = MesonOperatorEnv::registerAll();
  cout << "registered = " << registered << endl;

  XMLReader xml_in(argv[1]);
  FileTestParam params(xml_in, "/FileTest");

  XMLFileWriter xml_out(argv[2]);
  push(xml_out, "FileTest");

  Real   mass = 0.5234567;  // something non-trivial

  // Map directions onto names
  Array<string> txyz(4);
  txyz[0] = 't';
  txyz[1] = 'x';
  txyz[2] = 'y';
  txyz[3] = 'z';

  // Construct the source and sink suffices from the smearing state
  string source_suffix;
  if (params.source_smear[0] == 'P')
    source_suffix = "P";
  else if (params.source_smear[0] == 'D')
    source_suffix = "S";

  string sink_suffix;
  if (params.sink_smear[0] == 'P')
    sink_suffix = "P";
  else if (params.sink_smear[0] == 'D')
    sink_suffix = "S";

  // Grrh, map arbitrary names to arbitrary names
  // This is happening since the source-particle, sink-particle, etc.
  // are free parameters in chroma's mesonspec
  map<string,string> trans_map;
  trans_map["a0"] = "gamma0";
  trans_map["a0_x"] = "gamma0";
  trans_map["a0_y"] = "gamma0";
  trans_map["a0_z"] = "gamma0";

  trans_map["a1_x"] = "gamma14";
  trans_map["a1_y"] = "gamma13";
  trans_map["a1_z"] = "gamma11";

  trans_map["a2_x"] = "a2_x";
  trans_map["a2_y"] = "a2_y";
  trans_map["a2_z"] = "a2_z";

  trans_map["rho1_x"] = "gamma1";
  trans_map["rho1_y"] = "gamma2";
  trans_map["rho1_z"] = "gamma4";

  trans_map["b0"] = "b0";
  trans_map["b0_x"] = "b0_x";
  trans_map["b0_y"] = "b0_y";
  trans_map["b0_z"] = "b0_z";

  trans_map["b1_x"] = "gamma6";
  trans_map["b1_y"] = "gamma5";
  trans_map["b1_z"] = "gamma3";

  trans_map["b2_x"] = "b2_x";
  trans_map["b2_y"] = "b2_y";
  trans_map["b2_z"] = "b2_z";


 {

#if 0
   Array<int> p(Nd -1);
   p[0] = 1;
   p[1] = 0;
   p[2] = 0;
   cout << "polarisation vector" << endl;
    for(int dir_src=0; dir_src <= 3; ++dir_src)
      {
	for(int dir_snk=0; dir_snk <= 3; ++dir_snk)
	  {
	    cout << "dir_src = " << dir_src << endl;
	    cout << "dir_snk = " << dir_snk << endl;


	    // Do the polarization inner product
	    Complex sum = zero;
	    for(int r=0; r < 3; ++r)
	      {
		Complex Z_src =  (minkPolVec(mass, p, r, 36))[dir_src];
		Complex Z_snk =  (minkPolVec(mass, p, r, 36))[dir_snk];
		
		sum += conj(Z_snk)*Z_src;
		//cout << "term = " <<  conj(Z_snk)*Z_src << endl;
	      }

	    cout << "sqrt(mod2(sum)) = " << sqrt(mod2(sum)) << endl;
	  }
      }

    cout << "polarisation tensor" << endl;
    for(int dir_src=0; dir_src <= 3; ++dir_src)
      {
	for(int dir_snk=0; dir_snk <= 3; ++dir_snk)
	  {
	    cout << "dir_src = " << dir_src << endl;
	    cout << "dir_snk = " << dir_snk << endl;


	    // Do the polarization inner product
	    Complex sum = zero;
	    for(int r=0; r < 5; ++r)
	      {
		Complex Z_src =  (minkPolTens(mass, p, r, 36))[0][dir_src];
		Complex Z_snk =  (minkPolTens(mass, p, r, 36))[0][dir_snk];
		
		sum += conj(Z_snk)*Z_src;
		// cout << "term = " <<  conj(Z_snk)*Z_src << endl;
	      }

	    cout << "sqrt(mod2(sum)) = " << sqrt(mod2(sum)) << endl;
	  }
      }



#endif
 }



  vector<string> operator_list;
  operator_list.push_back("PIONxNABLA_T1");
  operator_list.push_back("A0xNABLA_T1");
  operator_list.push_back("A0_2xNABLA_T1");
  operator_list.push_back("RHOxNABLA_A1");
  operator_list.push_back("RHOxNABLA_T1");
  operator_list.push_back("RHOxNABLA_T2");
  operator_list.push_back("RHOxNABLA_E");
  // operator_list.push_back("A1xNABLA_A1");
  //operator_list.push_back("A1xNABLA_T1");
  operator_list.push_back("A1xNABLA_T2");
  operator_list.push_back("B1xNABLA_T1");
  operator_list.push_back("B1xNABLA_T2");
  //operator_list.push_back("A1xNABLA_E");
  //operator_list.push_back("A0_2xD_T2");
  //operator_list.push_back("RHOxD_A2");
  //operator_list.push_back("RHOxD_T1");
  //operator_list.push_back("RHOxD_T2");
  operator_list.push_back("A1xD_A2");
  operator_list.push_back("A1xD_T1");
  operator_list.push_back("A1xD_T2");
  operator_list.push_back("A1xD_E");
  //operator_list.push_back("RHOxB_T1");

  try
  {
    for(vector<string>::const_iterator op=operator_list.begin(); 
	op != operator_list.end(); 
	++op)
    {
      string source_id = *op;
      cout << "Source operator = " << source_id << endl;

      Handle< Operator<Complex> >  oper();

      //
      // Create the source operator
      //
      const string operator_path = "/FileTest";
	
      Handle< Operator<Complex> >
	sourceOperator(TheMesonOperatorFactory::Instance().createObject(source_id,
									xml_in,
									operator_path));

      string source_wavefunc = sourceOperator->operatorName();
      cout << "source_wavefunc = " << source_wavefunc << endl;

      string filename(source_wavefunc + string(".results"));
      ofstream of(filename.c_str());
   
      // Loop over source overlaps
      vector< Handle< WaveFunction<Complex> > > source_wvfs(sourceOperator->overlaps());

      for(vector< Handle< WaveFunction<Complex> > >::const_iterator src_wvf=source_wvfs.begin();
	  src_wvf != source_wvfs.end(); 
	  ++src_wvf)
      {
	const WaveFunction<Complex>& source_wvf = *(*src_wvf);
	string source_particle = source_wvf.particleName();

	cout << "source_particle = " << source_particle << endl;

	// Loop over simple Gamma(n) sink operators
	for(int n=0; n < 16; ++n)
	{
	  string sink_id;
	  {
	    ostringstream os;
	    os << "Gamma" << n;
	    sink_id = os.str();
	  }

	  //
	  // Create the sink operator
	  //
	  Handle< Operator<Complex> >
	    sinkOperator(TheMesonOperatorFactory::Instance().createObject(sink_id,
									  xml_in,
									  operator_path));

	  cout << "Sink operator = " << sink_id << endl;

	  // Loop over sink overlaps
	  vector< Handle< WaveFunction<Complex> > > sink_wvfs(sinkOperator->overlaps());

	  for(vector< Handle< WaveFunction<Complex> > >::const_iterator snk_wvf=sink_wvfs.begin();
	      snk_wvf != sink_wvfs.end(); 
	      ++snk_wvf)
	  {
	    const WaveFunction<Complex>& sink_wvf = *(*snk_wvf);
	    string sink_particle = sink_wvf.particleName();
	    string sink_wavefunc = sink_particle;

	    cout << "sink_wavefunc = " << sink_wavefunc << endl;
	    cout << "sink_particle = " << sink_particle << endl;
    
	    // Only if the quantum numbers match will we consider this cross-corr
	    if (source_particle != sink_particle)
	      continue;

	    cout << "Source and sink particle states match" << endl;
    
	    int dir_src_start=0, dir_src_end=0;
	    int dir_snk_start=0, dir_snk_end=0;

	    if (source_wvf.numDir() > 0)
	    {
	      dir_src_start = 1;
	      dir_src_end = source_wvf.numDir();
	    }

	    if (sink_wvf.numDir() > 0)
	    {
	      dir_snk_start = 1;
	      dir_snk_end = sink_wvf.numDir();
	    }

	    if (source_wvf.numPolar() != sink_wvf.numPolar())
	    {
	      cerr << argv[0] << ": something wrong, unequal number of source/sink polarizations" 
		   << endl;
	      exit(1);
	    }

	    int num_polar = max(source_wvf.numPolar(),1);   // source and sink the same quantum numbers

	    // Loop over all the momenta
	    int max_mom = params.mom2_max;
	    ArrayInt p(3);
	    for(p[0]=-max_mom; p[0] < max_mom; ++p[0])
	      for(p[1]=-max_mom; p[1] < max_mom; ++p[1])
		for(p[2]=-max_mom; p[2] < max_mom; ++p[2])
		{
		  int p_sq = norm2(p);
		  if (p_sq > params.mom2_max) continue;

		  printf("momentum = [%d,%d,%d]\n", p[0], p[1], p[2]);
    
		  // Loop over all the directions
		  for(int dir_src=dir_src_start; dir_src <= dir_src_end; ++dir_src)
		  {
		    for(int dir_snk=dir_snk_start; dir_snk <= dir_snk_end; ++dir_snk)
		    {
		      cout << "dir_src = " << dir_src << endl;
		      cout << "dir_snk = " << dir_snk << endl;
		      cout << "num_polar = " << num_polar << endl;

		      // Do the polarization inner product
		      Array<Complex> Z_src = source_wvf(mass, p, dir_src, 0);
		      Array<Complex> Z_snk = sink_wvf(mass, p, dir_snk, 0);
		      int len_src = Z_src.size();
		      int len_snk = Z_snk.size();
		      cout << "len_src = " << len_src << endl;		      
		      cout << "len_snk = " << len_snk << endl;		      

		      Array< Array<Complex> > Zsq(len_src);
		      for(int i=0; i < len_src; ++i)
			Zsq[i].resize(len_snk);

		      for(int i=0; i < len_src; ++i)
			for(int j=0; j < len_snk; ++j)
			  Zsq[i][j] = zero;
		       
		      // Complex sum = zero;
		      int rr = (num_polar - 1)/2; 
		      for(int r= - rr ; r < rr + 1; ++r)
		      {
			Array<Complex> Z_src = source_wvf(mass, p, dir_src, r);
			Array<Complex> Z_snk = sink_wvf(mass, p, dir_snk, r);
	      
			for(int i=0; i < len_src; ++i)
			  for(int j=0; j < len_snk; ++j)	
			    Zsq[i][j] += conj(Z_snk[j])*Z_src[i];
			
			//cout << "term = " <<  conj(Z_snk)*Z_src << endl;
		      }

		      //cout << "mod2(sum) = " << mod2(sum) << endl;
		      bool test = false;
		      for(int i=0; i < len_src; ++i)
			for(int j=0; j < len_snk; ++j)
			  {
			    Real Zsqsq = mod2( Zsq[i][j] );
			    cout << "Zsqsq = " <<  Zsqsq << endl;  
			    if  (toFloat( Zsqsq ) > params.tol)
			      {test = true;}
			  }

		      if (test)
		      {
			  cout << "writing out a filename " <<  endl; 
			// Construct the quantum number part of the file name
			string source_state = source_particle;
			if (dir_src > 0)
			{
			  source_state += "_" + txyz[dir_src];
			}

#if 0
			string sink_state = sink_particle;
			if (dir_snk > 0)
			{
			  sink_state += "_" + txyz[dir_snk];
			}
#else
			// HACK for now use the gamma name as the particle id
			string sink_state = sink_id;
#endif

			// If the source and sink quantum state are the same, then coalesce
			string meson_particle;
			if (source_state == sink_state)
			{
			  meson_particle = source_state;
			}
			else
			{
			  meson_particle = sink_state + "-" + source_state;
			}

			// Construct the total file name
			string name = meson_particle + mom_name(p) 
			  + "." + params.mass_string 
			  + "." + params.source_smear + "_" + "PT" 
			  + "." + params.sink_smear + "_" + source_wavefunc 
			  + "." + source_suffix + sink_suffix;
	
			//HARDWIRED FOR GAMMA MATRICES
			//string name = meson_particle + mom_name(p) 
			//  + "." + params.mass_string 
			//  + "." + params.source_smear + "_" + source_wavefunc 
			//  + "." + params.sink_smear + "_PT" 
			//  + "." + source_suffix + sink_suffix;
	    
			cout << name << endl;
			of << name << endl;
		      }
		    }
		  }
		}
	  }
	}
      }
    }
  }
  catch(std::bad_cast) 
  {
    cerr << argv[0] << ": caught cast error" << endl;
    exit(1);
  }
  catch(const std::string& e) 
  {
    cerr << argv[0] << ": Caught Exception: " << e << endl;
    exit(1);
  }
  catch(std::exception& e) 
  {
    cerr << argv[0] << ": Caught standard library exception: " << e.what() << endl;
    exit(1);
  }
  catch(...)
  {
    cerr << argv[0] << ": caught generic exception during measurement" << endl;
    exit(1);
  }

  pop(xml_out);  // FileTest
  xml_out.close();

  return 0;
}
