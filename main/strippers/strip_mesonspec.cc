// $Id: strip_mesonspec.cc,v 2.0 2008/12/05 04:44:04 edwards Exp $

#include "strippers.h"
#include <map>

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

  int          mom2_max;
  bool         avg_equiv_mom;

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


//! Names of sources
struct SourceSinkType_t
{
  string        source_type_1;
  string        sink_type_1;

  string        source_type_2;
  string        sink_type_2;
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
  Real             Mass;       // quark mass (bare units)
};


//! Propagator sink header
struct PropSink_t
{
  string           sink_type;   // Point, Shell, Wall, etc.
  SmearingParam_t  sinkSmearParam;
  // wvf-function smearing type (Gaussian, Exponential, etc.)
  // smearing width
  // number of iteration for smearing
  int              j_decay;         // Decay direction
};


struct Meson_momenta_t
{
  int sink_mom_num;
  Array<int> sink_mom;
  Array<Complex> mesprop;
};

struct Wilson_mesons_t
{
  Array<Meson_momenta_t> mom;
};

struct ForwardPropHeaders_t
{
  PropSource_t          source_header_1;
  PropSource_t          source_header_2;
  PropSink_t            sink_header_1;
  PropSink_t            sink_header_2;
};

struct Wilson_hadron_measurements_t
{
  int loop;

  string                       source_particle;
  string                       source_wavetype;
  string                       sink_particle;
  string                       sink_wavetype;

  Array<ForwardPropHeaders_t>  forward_prop_headers;
  Array<SourceSinkType_t>      source_sink_type;
  Real                         Mass_1;       // quark mass (bare units)
  Real                         Mass_2;       // quark mass (bare units)

  bool         Pt_src;   // Not read, but determined from source
  bool         Sl_src;   // Not read, but determined from source
  bool         Wl_src;   // Not read, but determined from source

  bool         Pt_snk;   // Not read, but determined from sink
  bool         Sl_snk;   // Not read, but determined from sink
  bool         Wl_snk;   // Not read, but determined from sink

  Wilson_mesons_t  meson;
};

struct MesonSpec_t
{
  Output_version_t output_version;
  Param_t          param;
  Array<Wilson_hadron_measurements_t>  had;
};


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
  
  ostringstream error_stream;

  int version;
  read(top, "version", version);

  switch (version)
  {
  case 1:
    break;

  default:
    error_stream << "Within Param read - unsupported param version = "
		 << version << endl;
    throw error_stream.str();
  }

  read(top, "mom2_max", param.mom2_max);
  read(top, "avg_equiv_mom", param.avg_equiv_mom);

  if(top.count("nrow")!=0)
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

//! Names of sources 
void read(XMLReader& xml, const string& path, SourceSinkType_t& param)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "source_type_1", param.source_type_1);
  read(paramtop, "sink_type_1", param.sink_type_1);

  read(paramtop, "source_type_2", param.source_type_2);
  read(paramtop, "sink_type_2", param.sink_type_2);
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
  break;

  default:
  {
    XMLReader sourcetop(paramtop, "Source");
    
    read(sourcetop, "SourceType", header.source_type);
    read(sourcetop, "j_decay",  header.j_decay);

    if (header.source_type == "SHELL_SOURCE" || header.source_type == "SF_SHELL_SOURCE" )
    {
      read(sourcetop, "SmearingParam", header.sourceSmearParam);
    }		
  }
  }
}


// Sink header read
void read(XMLReader& xml, const string& path, PropSink_t& header)
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
  {
    read(paramtop, "sink_type", header.sink_type);
    read(paramtop, "j_decay",  header.j_decay);

    if (header.sink_type == "SHELL_SINK")
    {
      XMLReader shelltop(paramtop, "ShellSink");
      read(shelltop, "SinkSmearingParam", header.sinkSmearParam);
    }
  }
  break;

  default:
  {
    XMLReader sinktop(paramtop, "Sink");
    
    read(sinktop, "SinkType", header.sink_type);
    read(sinktop, "j_decay",  header.j_decay);

    if (header.sink_type == "SHELL_SINK" || header.sink_type == "SF_SHELL_SINK")
    {
      read(sinktop, "SmearingParam", header.sinkSmearParam);
    }		
  }
  }
}


// Forward propagator header read
void read(XMLReader& xml, const string& path, ChromaProp_t& param)
{
  XMLReader paramtop(xml, path);
  int version;

  ostringstream error_stream;

  read(paramtop, "version", version);
  switch(version) {
  case 4: 
    read(paramtop, "Mass", param.Mass);
    break;
  case 5:
  case 6:
  case 7:
  case 8:
    read(paramtop, "FermionAction/Mass", param.Mass);
    break;
  default:
    error_stream << "Prop header version " << version << " unknown" << endl;
    throw error_stream.str();
  } 
}


// Read a momenta struct
void read(XMLReader& xml, const string& path, Meson_momenta_t& mom)
{
  XMLReader top(xml, path);

  read(top, "sink_mom_num", mom.sink_mom_num);
  read(top, "sink_mom", mom.sink_mom);
  read(top, "mesprop", mom.mesprop);
}


// meson struct
void read(XMLReader& xml, const string& path, Wilson_mesons_t& mes)
{
  XMLReader top(xml, path);

  read(top, "momenta", mes.mom);
}


// Source header read
void read(XMLReader& xml, const string& path, ForwardPropHeaders_t& header)
{
  XMLReader top(xml, path);

  read(top, "First_forward_prop/PropSource", header.source_header_1);
  if (top.count("First_forward_prop/PropSinkSmear") > 0)
    read(top, "First_forward_prop/PropSinkSmear", header.sink_header_1);
  else if (top.count("First_forward_prop/PropSink") > 0)
    read(top, "First_forward_prop/PropSink", header.sink_header_1);
  
  read(top, "Second_forward_prop/PropSource", header.source_header_2);
  if (top.count("Second_forward_prop/PropSinkSmear") > 0)
    read(top, "Second_forward_prop/PropSinkSmear", header.sink_header_2);
  else if (top.count("Second_forward_prop/PropSink") > 0)
    read(top, "Second_forward_prop/PropSink", header.sink_header_2);
}


// Read a hadron measurements
void read(XMLReader& xml, const string& path, Wilson_hadron_measurements_t& had)
{
//  cout << "Wilson_had = " << path << endl;

  XMLReader top(xml, path);

  read(top, "source_particle", had.source_particle);
  read(top, "source_wavetype", had.source_wavetype);
  read(top, "sink_particle", had.sink_particle);
  read(top, "sink_wavetype", had.sink_wavetype);

  read(top, "Mass_1", had.Mass_1);
  read(top, "Mass_2", had.Mass_2);

  read(top, "Forward_prop_headers", had.forward_prop_headers);
  read(top, "SourceSinkType", had.source_sink_type);

  // Determine what kind of source to use 
  had.Pt_src = (had.source_sink_type[0].source_type_1 == "POINT_SOURCE" ||
		had.source_sink_type[0].source_type_1 == "SF_POINT_SOURCE" ) ? true : false;

  had.Sl_src = (had.source_sink_type[0].source_type_1 == "SHELL_SOURCE" ||
		had.source_sink_type[0].source_type_1 == "SF_SHELL_SOURCE" ) ? true : false;

  had.Wl_src = (had.source_sink_type[0].source_type_1 == "WALL_SOURCE" ||
		had.source_sink_type[0].source_type_1 == "SF_WALL_SOURCE" ) ? true : false;

  // Determine what kind of sink to use 
  had.Pt_snk = (had.source_sink_type[0].sink_type_1 == "POINT_SINK" ||
		had.source_sink_type[0].sink_type_1 == "SF_POINT_SINK" ) ? true : false;

  had.Sl_snk = (had.source_sink_type[0].sink_type_1 == "SHELL_SINK" ||
		had.source_sink_type[0].sink_type_1 == "SF_SHELL_SINK" ) ? true : false;

  had.Wl_snk = (had.source_sink_type[0].sink_type_1 == "WALL_SINK" ||
		had.source_sink_type[0].sink_type_1 == "SF_WALL_SINK" ) ? true : false;

  string xpath;

  // Mesons
  xpath = "Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.meson);
}

// Mother of all readers
void read(XMLReader& xml, const string& path, MesonSpec_t& spec)
{
  ostringstream error_stream;

  try 
  {
    XMLReader mother(xml, path);

    read(mother, "Output_version", spec.output_version);

    read(mother, "Output_version", spec.output_version);
    if (mother.count("Input/spectrumOct_w/Param") != 0)
      read(mother, "Input/spectrumOct_w/Param", spec.param);
    else if (mother.count("Input/Param") != 0)
      read(mother, "Input/Param", spec.param);
    else
    {
      throw string("Cannot find Param input");
    }

    switch (spec.output_version.out_version)
    {
    case 2:
    case 3:
      break;

    default:
      error_stream << "read(mesonspec): Unsupported output version " 
		   << spec.output_version.out_version << endl;
      throw error_stream.str();
    }

    //Get the lattice size this way
    if(spec.param.nrow.size()==0)
      read(mother,"ProgramInfo/Setgeom/latt_size",spec.param.nrow);

    read(mother, "Wilson_hadron_measurements", spec.had);
  }
  catch( const string& error) { 
    cerr << "Error reading MesonSpec_t: " << error << endl;
    exit(1);
  }

}



/*
 * Similar structure to MesonSpec_t that holds filenames
 * NOTE: this cannot be merged into MesonSpec_t since that
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
  Array<File_t> mom;
};

struct File_hadron_measurements_t
{
  File_particles_t meson;
};

struct File_spectrum_t
{
  int nprop;
  Array<File_hadron_measurements_t>  had;
};



// Write baryons
ostream& operator<<(ostream& s, const Array<Complex>& p)
{
  for(int i=0; i < p.size(); ++i)
    s << i << " " << p[i].re << " " << p[i].im << "\n";
}


string mom_name(const Array<int>& sink_mom)
{
  ostringstream s;

  if (! (sink_mom[0] == 0 && sink_mom[1] == 0 && sink_mom[2] == 0) )
    s << "_px" << sink_mom[0] << "_py" << sink_mom[1] << "_pz" << sink_mom[2];

  return s.str();
}


// Write a complex array
void print_file(const string& filename, const Array<Complex>& barprop, 
		int nprop, 
		const MesonSpec_t& spec, 
		bool first)
{
  int j_decay = spec.had[0].forward_prop_headers[0].source_header_1.j_decay;

  if (first)
  {
    ofstream file(filename.c_str());

    // Write header line
    file << nprop << " " << spec.param.nrow[j_decay] << " 1 " 
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


void print_files(const File_spectrum_t& files, const MesonSpec_t& spec, bool first)
{
  cout << __func__ << ": entering" << endl;

  for(int i = 0; i < spec.had.size(); ++i) 
  {
    // mesons
    for(int n = 0; n < spec.had[i].meson.mom.size(); ++n)
      print_file(files.had[i].meson.mom[n].filename, 
		 spec.had[i].meson.mom[n].mesprop, 
		 files.nprop, spec, first);
  }

  cout << __func__ << ": exiting" << endl;
}


void construct_filenames(File_spectrum_t& files, const MesonSpec_t& spec,
			 int nprop)
{
  cout << __func__ << ": entering" << endl;

  files.nprop = nprop;   // the number of configurations

  // Resize to hold how many kappas are involved
  files.had.resize(spec.had.size());

  /*
   * Stringize names by looping once through array and finding
   * out things like momenta, etc.
   */

  int nontriv = -1;

  for(int i = 0; i < spec.had.size(); ++i) 
  {
    const Wilson_hadron_measurements_t& had = spec.had[i];

    // Sanity checks
    string source_meson_smear_state = had.source_wavetype;
    string sink_meson_smear_state = had.sink_wavetype;

    cout << "source meson state = " << had.source_particle << endl;
    cout << "sink meson state   = " << had.sink_particle << endl;

    // If the source and sink meson particle name is the same, make it just
    // one of the names. If not, then hypennate them.
    string meson_particle;
    if (had.source_particle == had.sink_particle)
    {
      meson_particle = had.source_particle;
    }
    else
    {
      meson_particle = had.source_particle + "-" + had.sink_particle;
    }

    cout << "meson_particle = " << meson_particle << endl;
    cout << "source_meson_smear_state = " << source_meson_smear_state << endl;
    cout << "sink_meson_smear_state = " << sink_meson_smear_state << endl;


    string Mass_1, Mass_2;
    string src_wvf_param_s;
    string snk_wvf_param_s;
    string source_string;
    string sink_string;
    string source_suffix, sink_suffix;

    {
      ostringstream Mass_str;
      Mass_str << int(10000*had.Mass_1 + 0.5);
      Mass_1 = Mass_str.str();
    }
    
    {
      ostringstream Mass_str;
      Mass_str << int(10000*had.Mass_2 + 0.5);
      Mass_2 = Mass_str.str();
    }
    
    /*
     * Source case
     */
    if (had.Pt_src)
    {
      source_string = ".P_" + source_meson_smear_state;
      source_suffix = "P";
    }
    else if (had.Wl_src)
    {
      source_string = ".W_" + source_meson_smear_state;
      source_suffix = "W";
    }
    else if (had.Sl_src)
    {
      // Width names. Do not want trailing zeros
      bool sameP = true;
      Real wvf_param = had.forward_prop_headers[0].source_header_1.sourceSmearParam.wvf_param;
      for(int m=0; m < had.forward_prop_headers.size(); ++m)
      {
	if (! ((had.forward_prop_headers[m].source_header_1.sourceSmearParam.wvf_param == wvf_param) &&
	       (had.forward_prop_headers[m].source_header_2.sourceSmearParam.wvf_param == wvf_param)))
	{
//	  cerr << "Not all source smearing widths are the same - cannot deal with that" << endl;
//	  exit(1);
	  sameP = false;
	}
      }

      if (sameP)
      {
	ostringstream wvf_param_str;

	wvf_param_str << "G";
	wvf_param_str << wvf_param;
	string wvf_ps = wvf_param_str.str();

	string::size_type idx = wvf_ps.find('.');  // replace '.' with 'p'
	if (idx != string::npos)
	  wvf_ps[idx] = 'p';

	src_wvf_param_s = wvf_ps;

	source_string = ".D" + src_wvf_param_s + "_" + source_meson_smear_state;
	source_suffix = "S";
      }
      else
      {
	// This is not handled in a clean way - hack for now
	ostringstream wvf_param_str_1;
	wvf_param_str_1 << "G";
	wvf_param_str_1 << had.forward_prop_headers[0].source_header_1.sourceSmearParam.wvf_param;
	string wvf_ps_1 = wvf_param_str_1.str();

	string::size_type idx_1 = wvf_ps_1.find('.');  // replace '.' with 'p'
	if (idx_1 != string::npos)
	  wvf_ps_1[idx_1] = 'p';

	ostringstream wvf_param_str_2;
	wvf_param_str_2 << "G";
	wvf_param_str_2 << had.forward_prop_headers[0].source_header_1.sourceSmearParam.wvf_param;
	string wvf_ps_2 = wvf_param_str_2.str();

	string::size_type idx_2 = wvf_ps_2.find('.');  // replace '.' with 'p'
	if (idx_2 != string::npos)
	  wvf_ps_2[idx_2] = 'p';

	source_string = ".H" + wvf_ps_1 + "-" + wvf_ps_2 + "_" + source_meson_smear_state;
	source_suffix = "S";
      }
    }
    else
    {
      cerr << __func__ << ": Unsupported source type" << endl;
      exit(1);
    }

    cout << "source_string = " << source_string << endl;
      
    /*
     * Sink case
     */
    if (had.Pt_snk)
    {
      sink_string = ".P_" + sink_meson_smear_state;
      sink_suffix = "P";
    }
    else if (had.Wl_snk)
    {
      sink_string = ".W_" + sink_meson_smear_state;
      sink_suffix = "W";
    }
    else if (had.Sl_snk)
    {
      // Width names. Do not want trailing zeros
      ostringstream wvf_param_str;

      Real wvf_param = had.forward_prop_headers[0].sink_header_1.sinkSmearParam.wvf_param;
      for(int m=0; m < had.forward_prop_headers.size(); ++m)
      {
	if (! ((had.forward_prop_headers[m].sink_header_1.sinkSmearParam.wvf_param == wvf_param) &&
	       (had.forward_prop_headers[m].sink_header_2.sinkSmearParam.wvf_param == wvf_param)))
	{
//	  cerr << "Not all sink smearing widths are the same - cannot deal with that" << endl;
//	  exit(1);
	}
      }

      wvf_param_str << "G";
      wvf_param_str << wvf_param;
      string wvf_ps = wvf_param_str.str();

      string::size_type idx = wvf_ps.find('.');  // replace '.' with 'p'
      if (idx != string::npos)
	wvf_ps[idx] = 'p';

      snk_wvf_param_s = wvf_ps;

      sink_string = ".D" + snk_wvf_param_s + "_" + sink_meson_smear_state;
      sink_suffix = "S";
    }
    else
    {
      cerr << __func__ << ": Unsupported sink type" << endl;
      exit(1);
    }

    cout << "sink_string = " << sink_string << endl;

    cout << "  Mass_1 = " << Mass_1 << endl
	 << "  Mass_2 = " << Mass_2 << endl;

    cout << "source_string = " << source_string << endl;
    cout << "sink_string = " << sink_string << endl;
    

    //
    // Construct meson names
    //
    cout << "Creating meson file names" << endl;

    {
      string Mass_12;
      if (Mass_1 == Mass_2)
      {
	Mass_12 = ".D" + Mass_1;
      }
      else
      {
	Mass_12 = ".H" + Mass_1+"_"+Mass_2;
      }

      files.had[i].meson.mom.resize(had.meson.mom.size());

      for(int n = 0; n < had.meson.mom.size(); ++n)
      {
	files.had[i].meson.mom[n].filename = meson_particle 
	  + mom_name(had.meson.mom[n].sink_mom) 
	  + Mass_12 + source_string + sink_string + "." + source_suffix + sink_suffix;
      }
    }
  }   // end for i

  cout << __func__ << ": exiting" << endl;
}


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

  string xml_in_root = "/MesonSpectrum";


  // Process the first file specially - read alls it params
  cout << "Open file " << argv[1] << endl;
  XMLReader xml_in(argv[1]);
  cout << "Done Open file " << argv[1] << endl;
  
  // Big nested structure that is image of entire file
  MesonSpec_t  spec;

  // Read data
  cout << "Read config 1 data " << argv[1] << endl;
  read(xml_in, xml_in_root, spec);
  xml_in.close();
  
  cout << "spec.had.size = " << spec.had.size() << endl;
  

  // We care about the output version
  switch (spec.output_version.out_version)
  {
  case 2:
  case 3:
    cout << "Processing output version " << spec.output_version.out_version << endl;
    break;

  default:
    cerr << "Unexpected output version " << spec.output_version.out_version << endl;
    exit(1);
  }


  // Structure holding file names
  File_spectrum_t files;

  // Construct the filenames into the structure  files
  cout << "Construct filenames" << endl;
  try {
    construct_filenames(files, spec, argc-1);
  }
  catch( const string& error) { 
    cerr << "Error constructing filenames: " << error << endl;
    exit(1);
  }

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

  cout << argv[0] << ": finished" << endl;

  exit(0);
}
