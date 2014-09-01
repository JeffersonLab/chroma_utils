// $Id: strip_spectrumQll.cc,v 2.0 2008/12/05 04:44:04 edwards Exp $
// WALL source or sink does not work

#include <iostream>
#include <fstream>

#include <cstdio>

#include "strippers.h"


using namespace Strippers;
using namespace std;



const int Ns = 4;
const int Nd = 4;

/*
 * Input 
 */
struct Output_version_t
{
  int out_version;
};

struct Lattis_t
{
};

struct Param_t
{
  int          version;

  bool         Pt_snk;
  bool         Sl_snk;
  bool         Wl_snk;
  string       wvf_kind;
  Array<Real>  wvf_param;
  Array<int>   wvfIntPar;

  Array<int>   nrow;
};

/*
 * Structures for hadron parts
 */
//! Parameters for sources and sinks
struct SmearingParam_t
{
  string        wvf_kind;
  Real          wvf_param;
  int           wvfIntPar;
};


// Relevant source header params
struct PropSource_t
{
  string           source_type;   // Point, Shell, Wall, etc.
  SmearingParam_t  sourceSmearParam;
                             // wvf-function smearing type (Gaussian, Exponential, etc.)
                             // smearing width
                             // number of iteration for smearing
  int              j_decay;         // Decay direction
};

// Relevant propagator inversion parameters
struct ChromaProp_t
{
  Real            Mass;       // quark mass (bare units)
};



struct Wilson_Qll_baryons_t
{
  Array<Complex> LambdaQ ;
  Array<Complex> SigmaQx ;
  Array<Complex> SigmaQy ;
  Array<Complex> SigmaQz ;
};

struct Wilson_Ql_mesons_t
{
  Array<Complex> HlPSmes ;
};

struct Wilson_hadron_measurements_t
{
  int loop;

  PropSource_t            source_header;
  ChromaProp_t            prop_header;

  bool         Pt_src;   // Not read, but determined from source
  bool         Sl_src;   // Not read, but determined from source
  bool         Wl_src;   // Not read, but determined from source

  Wilson_Qll_baryons_t  Qll_PP;
  Wilson_Qll_baryons_t  Qll_PS;
  Wilson_Qll_baryons_t  Qll_PW;
  

  Wilson_Qll_baryons_t  Qll_SP;
  Wilson_Qll_baryons_t  Qll_SS;
  Wilson_Qll_baryons_t  Qll_SW;

  Wilson_Ql_mesons_t   Ql_PP;
  Wilson_Ql_mesons_t   Ql_PS;
  Wilson_Ql_mesons_t   Ql_PW;

  Wilson_Ql_mesons_t   Ql_SP;
  Wilson_Ql_mesons_t   Ql_SS;
  Wilson_Ql_mesons_t   Ql_SW;
 

};

struct Spectrum_w_t
{
  Output_version_t output_version;
  Param_t          param;
  Array<Wilson_hadron_measurements_t>  had;
};


// Read a complex
void read(XMLReader& xml, const string& path, Complex& com)
{
  try {
    XMLReader complextop(xml, path);
    read(complextop, "re", com.re);
    read(complextop, "im", com.im);

  }
  catch( const string& error) { 
    cerr << "Error reading data: " << error << endl;
  }
}

/*
 * Read the parameters
 */
void read(XMLReader& xml, const string& path, Output_version_t& out)
{
  read(xml, path + "/out_version", out.out_version) ;
}

// params
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader top(xml, path);
  
  int version;
  read(top, "version", version);

  read(top, "Pt_snk", param.Pt_snk);
  read(top, "Sl_snk", param.Sl_snk);
 
  read(top, "wvf_kind", param.wvf_kind);
  read(top, "wvf_param", param.wvf_param);
  read(top, "wvfIntPar", param.wvfIntPar);

  if (param.wvf_param.size() != param.wvfIntPar.size())
  {
    cerr << "Unexpected wvf_param array size";
    throw;
  }

  read(top, "nrow", param.nrow);  // lattice size
}


//! Read a smearing param struct
void read(XMLReader& xml, const string& path, SmearingParam_t& param)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "wvf_kind", param.wvf_kind);
  read(paramtop, "wvf_param", param.wvf_param);
  read(paramtop, "wvfIntPar", param.wvfIntPar);
}

// Source header read
void read(XMLReader& xml, const string& path, PropSource_t& header)
{
  XMLReader paramtop(xml, path);

  int version;
  read(paramtop, "version", version);
  switch(version) 
    {
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
      {
	read(paramtop, "source_type", header.source_type);
	read(paramtop, "j_decay",  header.j_decay);
	
	if (header.source_type == "SHELL_SOURCE")
	  {
	    XMLReader shelltop(paramtop, "ShellSource");
	    read(shelltop, "SourceSmearingParam", header.sourceSmearParam);
	  }
      }
      break ;
    case 6:
      {
	XMLReader sourcetop(paramtop, "Source");
	read(sourcetop, "SourceType", header.source_type);
	read(sourcetop, "j_decay",  header.j_decay);
	if (header.source_type == "SHELL_SOURCE")
	  {
	    read(sourcetop, "SmearingParam", header.sourceSmearParam);
	  }
      }
      break ;
    default:
      {
	cerr << "Source header version " << version << " unknown" << endl;
	throw;
      }
     
    }

}


// Forward propagator header read
void read(XMLReader& xml, const string& path, ChromaProp_t& param)
{
  XMLReader paramtop(xml, path);
  int version;

  read(paramtop, "version", version);
  switch(version) {
  case 4: 
    read(paramtop, "Mass", param.Mass);
    break;
  case 5:
  case 6:
  case 7:
  case 8:
  case 9:
  case 10:
    read(paramtop, "FermionAction/Mass", param.Mass);
    break;
  default:
    cerr << "Prop header version " << version << " unknown" << endl;
    throw;
  } 
}



//  Wilson_Qll_baryons_t
void read(XMLReader& xml, const string& path, Wilson_Qll_baryons_t& Qll)
{
  XMLReader top(xml, path);

  read(top, "LambdaQ", Qll.LambdaQ);
  read(top, "SigmaQx", Qll.SigmaQx);
  read(top, "SigmaQy", Qll.SigmaQy);
  read(top, "SigmaQz", Qll.SigmaQz);

}

//  Wilson_Qll_mesons_t
void read(XMLReader& xml, const string& path, Wilson_Ql_mesons_t& Ql)
{
  XMLReader top(xml, path);

  read(top, "HeavyLight", Ql.HlPSmes);
}

// Read a hadron measurements
void read(XMLReader& xml, const string& path, Wilson_hadron_measurements_t& had)
{
//  cout << "Wilson_had = " << path << endl;

  XMLReader top(xml, path);

  read(top, "loop", had.loop);

  read(top, "PropSource", had.source_header);
  read(top, "ForwardProp", had.prop_header);

  // Determine what kind of source to use
  had.Pt_src = (had.source_header.source_type == "POINT_SOURCE") ? true : false;
  had.Sl_src = (had.source_header.source_type == "SHELL_SOURCE") ? true : false;
  had.Wl_src = (had.source_header.source_type == "WALL_SOURCE")  ? true : false;

  string xpath;

  // 
  xpath = "Point_Point_Wilson_QllBaryons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.Qll_PP);

  xpath = "Point_Shell_Wilson_QllBaryons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.Qll_PS);

  xpath = "Point_Wall_Wilson_QllBaryons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.Qll_PW);

  // 
  xpath = "Shell_Point_Wilson_QllBaryons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.Qll_SP);

  xpath = "Shell_Shell_Wilson_QllBaryons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.Qll_SS);

  xpath = "Shell_Wall_Wilson_QllBaryons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.Qll_SW);

  xpath = "Point_Point_Wilson_Qlmeson";
  if (top.count(xpath) != 0)
    read(top, xpath, had.Ql_PP);
  
  xpath = "Point_Shell_Wilson_Qlmeson";
  if (top.count(xpath) != 0)
    read(top, xpath, had.Ql_PS);
  
  xpath = "Point_Wall_Wilson_Qlmeson";
  if (top.count(xpath) != 0)
    read(top, xpath, had.Ql_PW);
  
  xpath = "Shell_Point_Wilson_Qlmeson";
  if (top.count(xpath) != 0)
    read(top, xpath, had.Ql_SP);
  
  xpath = "Shell_Shell_Wilson_Qlmeson";
  if (top.count(xpath) != 0)
    read(top, xpath, had.Ql_SS);
  
  xpath = "Shell_Wall_Wilson_Qlmeson";
  if (top.count(xpath) != 0)
    read(top, xpath, had.Ql_SW);
  
}

// Mother of all readers
void read(XMLReader& xml, const string& path, Spectrum_w_t& spec)
{

  try 
  {
    XMLReader mother(xml, path);

    read(mother, "Output_version", spec.output_version);

    //

    switch (spec.output_version.out_version)
    {
    case 1:  
      //cout<<" Doing output version 1"<<endl ;
      read(mother, "Input/spectrumQll_w/Param", spec.param);
      read(mother, "Wilson_hadron_measurements", spec.had);
      break;
    case 2:
      //cout<<" Doing output version 2"<<endl ;
      read(mother, "Input/Param", spec.param);
      read(mother, "Wilson_hadron_measurements", spec.had);
      break;
    default:
      cerr << "read(spectrum): Unsupported output version " 
	   << spec.output_version.out_version<< endl;
      throw;
    }

  }
  catch( const string& error) { 
    cerr << "Error reading data: " << error << endl;
    exit(1);
  }

}



/*
 * Similar structure to Spectrum_w_t that holds filenames
 * NOTE: this cannot be merged into Spectrum_w_t since that
 * is volatile - upon XML reading this struct old data is
 * wiped out
 */
/*
 * Structures for hadron parts
 */
struct File_t
{
  string    filename;   // used by file writer below
};

struct File_particles_t
{
  File_t LambdaQ;
  File_t SigmaQx;
  File_t SigmaQy;
  File_t SigmaQz;
  File_t HlPSmes;
};


struct File_hadron_measurements_t
{
  File_particles_t Qll_PP;
  File_particles_t Qll_PS;
  File_particles_t Qll_PW;
  File_particles_t Qll_SP;
  File_particles_t Qll_SS;
  File_particles_t Qll_SW;

  File_particles_t Ql_PP ;
  File_particles_t Ql_PS ;
  File_particles_t Ql_PW ;
  File_particles_t Ql_SP ;
  File_particles_t Ql_SS ;
  File_particles_t Ql_SW ;

};

struct File_spectrum_t
{
  int nprop;
  Array<File_hadron_measurements_t>  had;
};



// Write mesons
ostream& operator<<(ostream& s, const Array<Real>& p)
{
  for(int i=0; i < p.size(); ++i)
    s << i << " " << p[i] << "\n";
}


// Write baryons
ostream& operator<<(ostream& s, const Array<Complex>& p)
{
  for(int i=0; i < p.size(); ++i)
    s << i << " " << p[i].re << "\n";
}


// Complex writer for now only writes real part   !!!!
void print_file(const string& filename, const Array<Complex>& barprop, 
		int nprop, 
		const Spectrum_w_t& spec, 
		bool first)
{
  int j_decay = spec.had[0].source_header.j_decay;

  //cout<<"Writting file: |"<<filename<<"|"<<endl ;
  if (first)
  {
    ofstream file(filename.c_str());

    // Write header line
    file << nprop << " " << spec.param.nrow[j_decay] << " 0 " 
	 << spec.param.nrow[0] << " 1" << endl;

    file << barprop;   // Write the file
    file.close();
  }
  else
  {
    ofstream file(filename.c_str(), std::ios_base::app);
    file << barprop;
    file.close();
  }
}


void print_files(const File_spectrum_t& files, const Spectrum_w_t& spec, bool first)
{
  for(int i = 0; i < spec.had.size(); ++i) 
    {
      if ( spec.had[i].Pt_src && spec.param.Pt_snk ) 
	{
	  print_file(files.had[i].Qll_PP.LambdaQ.filename, 
		     spec.had[i].Qll_PP.LambdaQ, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Qll_PP.SigmaQx.filename, 
		     spec.had[i].Qll_PP.SigmaQx, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Qll_PP.SigmaQy.filename, 
		     spec.had[i].Qll_PP.SigmaQy, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Qll_PP.SigmaQz.filename, 
		     spec.had[i].Qll_PP.SigmaQz, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Ql_PP.HlPSmes.filename, 
		     spec.had[i].Ql_PP.HlPSmes, 
		     files.nprop, spec, first);
	}
      
    
  
      if ( spec.had[i].Pt_src && spec.param.Sl_snk ) 
	{
	  print_file(files.had[i].Qll_PS.LambdaQ.filename, 
		     spec.had[i].Qll_PS.LambdaQ, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Qll_PS.SigmaQx.filename, 
		     spec.had[i].Qll_PS.SigmaQx, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Qll_PS.SigmaQy.filename, 
		     spec.had[i].Qll_PS.SigmaQy, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Qll_PS.SigmaQz.filename, 
		     spec.had[i].Qll_PS.SigmaQz, 
		     files.nprop, spec, first); 
	  print_file(files.had[i].Ql_PS.HlPSmes.filename, 
		     spec.had[i].Ql_PS.HlPSmes, 
		     files.nprop, spec, first);
	}
      
      if ( spec.had[i].Pt_src && spec.param.Wl_snk ) 
	{
	  print_file(files.had[i].Qll_PW.LambdaQ.filename, 
		     spec.had[i].Qll_PW.LambdaQ, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Qll_PW.SigmaQx.filename, 
		   spec.had[i].Qll_PW.SigmaQx, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Qll_PW.SigmaQy.filename, 
		     spec.had[i].Qll_PW.SigmaQy, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Qll_PW.SigmaQz.filename, 
		     spec.had[i].Qll_PW.SigmaQz, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Ql_PW.HlPSmes.filename, 
		     spec.had[i].Ql_PW.HlPSmes, 
		     files.nprop, spec, first);
	}

      if ( spec.had[i].Sl_src && spec.param.Pt_snk ) 
	{
	  print_file(files.had[i].Qll_SP.LambdaQ.filename, 
		     spec.had[i].Qll_SP.LambdaQ, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Qll_SP.SigmaQx.filename, 
		     spec.had[i].Qll_SP.SigmaQx, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Qll_SP.SigmaQy.filename, 
		     spec.had[i].Qll_SP.SigmaQy, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Qll_SP.SigmaQz.filename, 
		     spec.had[i].Qll_SP.SigmaQz, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Ql_SP.HlPSmes.filename, 
		     spec.had[i].Ql_SP.HlPSmes, 
		     files.nprop, spec, first);
	}

      

      if ( spec.had[i].Sl_src && spec.param.Sl_snk ) 
	{
	  print_file(files.had[i].Qll_SS.LambdaQ.filename, 
		     spec.had[i].Qll_SS.LambdaQ, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Qll_SS.SigmaQx.filename, 
		     spec.had[i].Qll_SS.SigmaQx, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Qll_SS.SigmaQy.filename, 
		     spec.had[i].Qll_SS.SigmaQy, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Qll_SS.SigmaQz.filename, 
		     spec.had[i].Qll_SS.SigmaQz, 
		     files.nprop, spec, first);
	  print_file(files.had[i].Ql_SS.HlPSmes.filename, 
		     spec.had[i].Ql_SS.HlPSmes, 
		     files.nprop, spec, first);
	} 
    
    if ( spec.had[i].Sl_src && spec.param.Wl_snk ) 
      {
	print_file(files.had[i].Qll_SW.LambdaQ.filename, 
		   spec.had[i].Qll_SW.LambdaQ, 
		   files.nprop, spec, first);
	print_file(files.had[i].Qll_SW.SigmaQx.filename, 
		   spec.had[i].Qll_SW.SigmaQx, 
		   files.nprop, spec, first);
	print_file(files.had[i].Qll_SW.SigmaQy.filename, 
		   spec.had[i].Qll_SW.SigmaQy, 
		   files.nprop, spec, first);
	print_file(files.had[i].Qll_SW.SigmaQz.filename, 
		   spec.had[i].Qll_SW.SigmaQz, 
		   files.nprop, spec, first);
	print_file(files.had[i].Ql_SW.HlPSmes.filename, 
		   spec.had[i].Ql_SW.HlPSmes, 
		   files.nprop, spec, first);
      } 
    }
}


void construct_filenames(File_spectrum_t& files, const Spectrum_w_t& spec,
			 int nprop)
{
  files.nprop = nprop;   // the number of configurations

  // Resize to hold how many kappas are involved
  files.had.resize(spec.had.size());

  /*
   * Stringize names by looping once through array and finding
   * out things like momenta, etc.
   */
  Array<string> Mass_s(spec.had.size());
  Array<string> src_wvf_param_s(spec.had.size());
  Array<string> snk_wvf_param_s(spec.had.size());

  for(int i = 0; i < spec.had.size(); ++i) 
  {
    // Mass values as an integer
    ostringstream Mass_str;
    Mass_str << int(10000*spec.had[i].prop_header.Mass + 0.5);
    Mass_s[i] = Mass_str.str();

    /*
     * Source case
     */
    if (spec.had[i].Sl_src)
    {
      // Width names. Do not want trailing zeros
      ostringstream wvf_param_str;

      if (spec.had[i].source_header.sourceSmearParam.wvf_kind == "GAUGE_INV_GAUSSIAN") {
	wvf_param_str << "G";
      } else {
	cerr << "Unknown or unsupported source wvf_kind" ;
	throw;
      }

      wvf_param_str << spec.had[i].source_header.sourceSmearParam.wvf_param;
      string wvf_ps = wvf_param_str.str();

      string::size_type idx = wvf_ps.find('.');  // replace '.' with 'p'
      if (idx != string::npos)
	wvf_ps[idx] = 'p';

      src_wvf_param_s[i] = wvf_ps;
    }
    else
    {
      src_wvf_param_s[i] = "";
    }
      
    
    /*
     * Sink case
     */
    if (spec.param.Sl_snk)
    {
      // Width names. Do not want trailing zeros
      ostringstream wvf_param_str;

      if (spec.param.wvf_kind == "") {
	wvf_param_str << "G";
      }
      else if (spec.param.wvf_kind == "GAUGE_INV_GAUSSIAN") {
	wvf_param_str << "G";
      } else {
	cerr << "Unknown or unsupported sink wvf_kind" ;
	throw;
      }

      wvf_param_str << spec.param.wvf_param[i];
      string wvf_ps = wvf_param_str.str();

      string::size_type idx = wvf_ps.find('.');  // replace '.' with 'p'
      if (idx != string::npos)
	wvf_ps[idx] = 'p';

      snk_wvf_param_s[i] = wvf_ps;
    }
    else
    {
      snk_wvf_param_s[i] = "";
    }
    
    
    cout << "  Mass_s[" << i << "]= " <<  Mass_s[i];
    if (spec.had[i].Sl_src)
      cout << "  src_wvf_param_s[" << i << "]= " << src_wvf_param_s[i];
    if (spec.param.Sl_snk)
      cout << "  snk_wvf_param_s[" << i << "]= " << snk_wvf_param_s[i];
    cout << endl;
  }  // end for (i)


    /*
     * Create source-sink smearing names
     */
    cout << "Creating file names" << endl;

    // For now assume only diagonal valence mesons
    for(int i = 0; i < spec.had.size(); ++i) 
    {
      /*
      cout<<"PT source: "<< spec.had[i].Pt_src<<endl ; 
      cout<<"SL source: "<< spec.had[i].Sl_src<<endl ;

      cout<<"PT sink: "<< spec.param.Pt_snk<<endl ;
      cout<<"SL sink: "<< spec.param.Sl_snk<<endl ;
      cout<<"WL sink: "<< spec.param.Wl_snk<<endl ;
      */

      int i1=i, i2=i;
      bool diag = true;
      string Mass_12 = ".D" + Mass_s[i1];
      string wvf_param12_s = "";   // ignore for the moment, should fill in

      //cout<<"CASE PP: "<< ( spec.had[i].Pt_src && spec.param.Pt_snk ) <<endl ;
      if ( spec.had[i].Pt_src && spec.param.Pt_snk ) 
	{
	  files.had[i].Qll_PP.LambdaQ.filename = "LambdaQ"+Mass_12+".P"+".P"+".PP" ;
	  files.had[i].Qll_PP.SigmaQx.filename = "SigmaQx"+Mass_12+".P"+".P"+".PP" ;
	  files.had[i].Qll_PP.SigmaQy.filename = "SigmaQy"+Mass_12+".P"+".P"+".PP" ;
	  files.had[i].Qll_PP.SigmaQz.filename = "SigmaQz"+Mass_12+".P"+".P"+".PP" ;
	  files.had[i].Ql_PP.HlPSmes.filename = "HlPSmes"+Mass_12+".P"+".P"+".PP" ;
	} 

      //cout<<"CASE PS: "<< ( spec.had[i].Pt_src && spec.param.Sl_snk )  <<endl ;
      if ( spec.had[i].Pt_src && spec.param.Sl_snk ) 
	{
	  files.had[i].Qll_PS.LambdaQ.filename = "LambdaQ"+Mass_12+".P"+".D"+ snk_wvf_param_s[i2]+".PS" ;
	  files.had[i].Qll_PS.SigmaQx.filename = "SigmaQx"+Mass_12+".P"+".D"+ snk_wvf_param_s[i2]+".PS" ;
	  files.had[i].Qll_PS.SigmaQy.filename = "SigmaQy"+Mass_12+".P"+".D"+ snk_wvf_param_s[i2]+".PS" ;
	  files.had[i].Qll_PS.SigmaQz.filename = "SigmaQz"+Mass_12+".P"+".D"+ snk_wvf_param_s[i2]+".PS" ;
	  files.had[i].Ql_PS.HlPSmes.filename = "HlPSmes"+Mass_12+".P"+".D"+ snk_wvf_param_s[i2]+".PS" ;
	}


      //cout<<"CASE SP: "<< ( spec.had[i].Sl_src && spec.param.Pt_snk )  <<endl ;
      if ( spec.had[i].Sl_src && spec.param.Pt_snk ) 
	{
	  files.had[i].Qll_SP.LambdaQ.filename = "LambdaQ"+Mass_12+".D"+ snk_wvf_param_s[i1]+".P"+".SP" ;
	  files.had[i].Qll_SP.SigmaQx.filename = "SigmaQx"+Mass_12+".D"+ snk_wvf_param_s[i1]+".P"+".SP" ;
	  files.had[i].Qll_SP.SigmaQy.filename = "SigmaQy"+Mass_12+".D"+ snk_wvf_param_s[i1]+".P"+".SP" ;
	  files.had[i].Qll_SP.SigmaQz.filename = "SigmaQz"+Mass_12+".D"+ snk_wvf_param_s[i1]+".P"+".SP" ;
	  files.had[i].Ql_SP.HlPSmes.filename = "HlPSmes"+Mass_12+".D"+ snk_wvf_param_s[i1]+".P"+".SP" ;
	}

      //cout<<"CASE SS: "<< ( spec.had[i].Sl_src && spec.param.Sl_snk )  <<endl ;
      if ( spec.had[i].Sl_src && spec.param.Sl_snk ) 
	{
	  files.had[i].Qll_SS.LambdaQ.filename = "LambdaQ"+Mass_12+".D"+ snk_wvf_param_s[i1]+".D"+ snk_wvf_param_s[i2]+".SS" ;
	  files.had[i].Qll_SS.SigmaQx.filename = "SigmaQx"+Mass_12+".D"+ snk_wvf_param_s[i1]+".D"+ snk_wvf_param_s[i2]+".SS" ;
	  files.had[i].Qll_SS.SigmaQy.filename = "SigmaQy"+Mass_12+".D"+ snk_wvf_param_s[i1]+".D"+ snk_wvf_param_s[i2]+".SS" ;
	  files.had[i].Qll_SS.SigmaQz.filename = "SigmaQz"+Mass_12+".D"+ snk_wvf_param_s[i1]+".D"+ snk_wvf_param_s[i2]+".SS" ;
	  files.had[i].Ql_SS.HlPSmes.filename = "HlPSmes"+Mass_12+".D"+ snk_wvf_param_s[i1]+".D"+ snk_wvf_param_s[i2]+".SS" ;
	}

      
    }  // end for i
    //DEBUG
    /**
    for(int i = 0; i < spec.had.size(); ++i) 
      {
	cout<<"FILE: "<<files.had[i].Qll_PP.LambdaQ.filename<<endl ;
	cout<<"FILE: "<<files.had[i].Qll_PS.LambdaQ.filename<<endl ;
	cout<<"FILE: "<<files.had[i].Qll_SS.LambdaQ.filename<<endl ;
	cout<<"FILE: "<<files.had[i].Qll_SP.LambdaQ.filename<<endl ;
      }
    **/
}

// I AM HERE !   !!!!!!!!!!!!!!!!!!!!!!!!!

//
// Main program - loop over files
//
int main(int argc, char **argv)
{
  if (argc == 1)
  {
    cerr << "Usage: " << argv[0] << " file1 [file2 file3... fileN]" << endl;
    exit(1);
  }

  string xml_in_root = "/spectrumQll_w";


  // Process the first file specially - read alls it params
  cout << "Open file " << argv[1] << endl;
  XMLReader xml_in(argv[1]);
  cout << "Done Open file " << argv[1] << endl;
  
  // Big nested structure that is image of entire file
  Spectrum_w_t  spec;

  // Read data
  cout << "Read config 1 data " << argv[1] << endl;
  read(xml_in, xml_in_root, spec);
  xml_in.close();
  
  cout << "spec.had.size = " << spec.had.size() << endl;
  cout << "spec.had.Qll_PP.Lambda.size = " << spec.had[0].Qll_PP.LambdaQ.size() << endl;
  cout << "spec.had.Qll_PS.Lambda.size = " << spec.had[0].Qll_PS.LambdaQ.size() << endl;
  cout << "spec.had.Qll_PW.Lambda.size = " << spec.had[0].Qll_PW.LambdaQ.size() << endl;

  cout << "spec.had.Qll_SP.Lambda.size = " << spec.had[0].Qll_SP.LambdaQ.size() << endl;
  cout << "spec.had.Qll_SS.Lambda.size = " << spec.had[0].Qll_SS.LambdaQ.size() << endl;
  cout << "spec.had.Qll_SW.Lambda.size = " << spec.had[0].Qll_SW.LambdaQ.size() << endl;
  
  if (spec.param.wvf_param.size() != spec.had.size())
  {
    cerr << "wvf_param and Mass array do not agree in size";
    throw;
  }


  // We care about the output version
  switch (spec.output_version.out_version)
  {
  case 1:
  case 2:
    cout << "Processing output version " << spec.output_version.out_version << endl;
    break;

  default:
    cerr << "Unexpected output version " << spec.output_version.out_version << endl;
    throw;
  }


  // Structure holding file names
  File_spectrum_t files;

  // Construct the filenames into the structure  files
  cout << "Construct filenames" << endl;
  construct_filenames(files, spec, argc-1);

  // Print the first time the files
  print_files(files, spec, true);

  for(int nfile=2; nfile < argc; ++nfile)
  {
    // Read data
    /* NOTE: could just read the hadron part */
    cout << "Reading config " << nfile << " - " << argv[nfile] << endl;
    xml_in.open(argv[nfile]);
    read(xml_in, xml_in_root, spec);
    xml_in.close();

    // Append to the files
    print_files(files, spec, false);
  }

  exit(0);
}
