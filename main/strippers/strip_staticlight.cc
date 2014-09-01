// $Id: strip_staticlight.cc,v 2.0 2008/12/05 04:44:04 edwards Exp $

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

  bool         MesonP;
  bool         BaryonP;

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
  Real            Mass;       // quark mass (bare units)
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



struct Wilson_mesons_t
{
  Array<Complex> mes1;
  Array<Complex> mes2; // when masses are different
};

struct Wilson_baryons_t
{
  Array<Complex> Lambda;
  Array<Complex> SigmaX;
  Array<Complex> SigmaY;
  Array<Complex> SigmaZ;
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

  ForwardPropHeaders_t  forward_prop_headers;
  SourceSinkType_t      source_sink_type;
  Real                  Mass_1;       // quark mass (bare units)
  Real                  Mass_2;       // quark mass (bare units)

  bool         Pt_src;   // Not read, but determined from source
  bool         Sl_src;   // Not read, but determined from source
  bool         Wl_src;   // Not read, but determined from source

  bool         Pt_snk;   // Not read, but determined from sink
  bool         Sl_snk;   // Not read, but determined from sink
  bool         Wl_snk;   // Not read, but determined from sink

  Wilson_mesons_t  meson;
  Wilson_baryons_t baryon;
};

struct HadSpec_t
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

  read(top, "MesonP", param.MesonP);
  read(top, "BaryonP", param.BaryonP);

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

    if (header.source_type == "SHELL_SOURCE")
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

    if (header.sink_type == "SHELL_SINK")
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


void read(XMLReader& xml, const string& path, Wilson_mesons_t& mes)
{
  XMLReader top(xml, path);

  read(top, "HeavyLight", mes.mes1);
  mes.mes2=mes.mes1 ;
  
}

// baryon struct
void read(XMLReader& xml, const string& path, Wilson_baryons_t& bar)
{
  XMLReader top(xml, path);

  read(top, "LambdaQ", bar.Lambda);
  read(top, "SigmaQx", bar.SigmaX);
  read(top, "SigmaQy", bar.SigmaY);
  read(top, "SigmaQz", bar.SigmaZ);
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

  read(top, "Mass_1", had.Mass_1);
  read(top, "Mass_2", had.Mass_2);

  read(top, "Forward_prop_headers", had.forward_prop_headers);
  read(top, "SourceSinkType", had.source_sink_type);

  string source_type;
  string sink_type;

  // Determine what kind of source to use 
  string pt_src = "POINT_SOURCE";
  if (had.source_sink_type.source_type_1 == pt_src)
  {
    if (! ((had.source_sink_type.source_type_1 == pt_src) &&
	   (had.source_sink_type.source_type_2 == pt_src)))
    {
      cerr << "Wilson_hadron_measurements_t: not all sources are a point_source" << endl;
      exit(1);
    }
    had.Pt_src = true;
    source_type = "Point";
  }
  else
  {
    had.Pt_src = false;
  }

  string sh_src = "SHELL_SOURCE";
  if (had.source_sink_type.source_type_1 == sh_src)
  {
    if (! ((had.source_sink_type.source_type_1 == sh_src) &&
	   (had.source_sink_type.source_type_2 == sh_src)))
    {
      cerr << "Wilson_hadron_measurements_t: not all sources are a shell_source" << endl;
      exit(1);
    }
    had.Sl_src = true;
    source_type = "Shell";
  }
  else
  {
    had.Sl_src = false;
  }

  string wl_src = "WALL_SOURCE";
  if (had.source_sink_type.source_type_1 == wl_src)
  {
    if (! ((had.source_sink_type.source_type_1 == wl_src) &&
	   (had.source_sink_type.source_type_2 == wl_src)))
    {
      cerr << "Wilson_hadron_measurements_t: not all sources are a wall_source" << endl;
      exit(1);
    }
    had.Wl_src = true;
    source_type = "Wall";
  }
  else
  {
    had.Wl_src = false;
  }


  // Determine what kind of sink to use 
  string pt_snk = "POINT_SINK";
  if (had.source_sink_type.sink_type_1 == pt_snk)
  {
    if (! ((had.source_sink_type.sink_type_1 == pt_snk) &&
	   (had.source_sink_type.sink_type_2 == pt_snk)))
    {
      cerr << "Wilson_hadron_measurements_t: not all sinks are a point_sink" << endl;
      exit(1);
    }
    had.Pt_snk = true;
    sink_type = "Point";
  }
  else
  {
    had.Pt_snk = false;
  }

  string sh_snk = "SHELL_SINK";
  if (had.source_sink_type.sink_type_1 == sh_snk)
  {
    if (! ((had.source_sink_type.sink_type_1 == sh_snk) &&
	   (had.source_sink_type.sink_type_2 == sh_snk)))
    {
      cerr << "Wilson_hadron_measurements_t: not all sinks are a shell_sink" << endl;
      exit(1);
    }
    had.Sl_snk = true;
    sink_type = "Shell";
  }
  else
  {
    had.Sl_snk = false;
  }

  string wl_snk = "WALL_SINK";
  if (had.source_sink_type.sink_type_1 == wl_snk)
  {
    if (! ((had.source_sink_type.sink_type_1 == wl_snk) &&
	   (had.source_sink_type.sink_type_2 == wl_snk)))
    {
      cerr << "Wilson_hadron_measurements_t: not all sinks are a wall_sink" << endl;
      exit(1);
    }
    had.Wl_snk = true;
    sink_type = "Wall";
  }
  else
  {
    had.Wl_snk = false;
  }

  string xpath;

  // Mesons
  xpath = source_type + "_" + sink_type + "_Wilson_Qlmeson";
  if (top.count(xpath) != 0){
    if(top.count(xpath) == 2){
      Wilson_mesons_t tt ;
      string xx = xpath+"[1]" ;
      read(top, xx, had.meson);
      xx = xpath+"[2]" ;
      read(top, xx, tt);
      had.meson.mes2 = tt.mes2 ;
    }
    else{
      read(top, xpath, had.meson);
    }
  }

    


  // Baryons
  xpath = source_type + "_" + sink_type + "_Wilson_QllBaryons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.baryon);
}

// Mother of all readers
void read(XMLReader& xml, const string& path, HadSpec_t& spec)
{
  ostringstream error_stream;

  try 
  {
    XMLReader mother(xml, path);

    read(mother, "Output_version", spec.output_version);

    read(mother, "Output_version", spec.output_version);
    if (mother.count("Input/Param") != 0)
      read(mother, "Input/Param", spec.param);
    else
    {
      throw string("Cannot find Param input");
    }

    switch (spec.output_version.out_version)
    {
    case 14:
    case 15:
      break;

    default:
      error_stream << "read(spectrum): Unsupported output version " 
		   << spec.output_version.out_version << endl;
      throw error_stream.str();
    }

    //Get the lattice size this way
    read(mother,"ProgramInfo/Setgeom/latt_size",spec.param.nrow);

    read(mother, "Wilson_hadron_measurements", spec.had);
  }
  catch( const string& error) { 
    cerr << "Error reading HadSpec_t: " << error << endl;
    exit(1);
  }

}



/*
 * Similar structure to HadSpec_t that holds filenames
 * NOTE: this cannot be merged into HadSpec_t since that
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

struct File_meson_t
{
  File_t    mes1;   // used by file writer below
  File_t    mes2;   // used by file writer below
};

struct File_baryon_t
{
  File_t    Lambda;   // used by file writer below
  File_t    SigmaX;   // used by file writer below
  File_t    SigmaY;   // used by file writer below
  File_t    SigmaZ;   // used by file writer below
};


struct File_hadron_measurements_t
{
  File_meson_t  meson;
  File_baryon_t baryon;
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
    s << i << " " << p[i].re << " " << p[i].im << "\n";
}


string mom_name(const Array<int>& sink_mom)
{
  ostringstream s;

  if (! (sink_mom[0] == 0 && sink_mom[1] == 0 && sink_mom[2] == 0) )
    s << "_px" << sink_mom[0] << "_py" << sink_mom[1] << "_pz" << sink_mom[2];

  return s.str();
}


void print_file(const string& filename, const Array<Real>& prop, 
		int nprop, 
		const HadSpec_t& spec, 
		bool first)
{
  int j_decay = spec.had[0].forward_prop_headers.source_header_1.j_decay;

  if (first)
  {
    ofstream file(filename.c_str());

    // Write header line
    file << nprop << " " << spec.param.nrow[j_decay] << " 0 " 
	 << spec.param.nrow[0] << " 1" << endl;

    file << prop;   // Write the file
    file.close();
  }
  else
  {
    ofstream file(filename.c_str(), std::ios_base::app);
    file << prop;
    file.close();
  }
}


// Complex writer
void print_file(const string& filename, const Array<Complex>& prop, 
		int nprop, 
		const HadSpec_t& spec, 
		bool first)
{
  int j_decay = spec.had[0].forward_prop_headers.source_header_1.j_decay;

  if (first)
  {
    ofstream file(filename.c_str());

    // Write header line
    file << nprop << " " << spec.param.nrow[j_decay] << " 1 " 
	 << spec.param.nrow[0] << " 1" << endl;

    file << prop;   // Write the file
    file.close();
  }
  else
  {
    ofstream file(filename.c_str(), std::ios_base::app);
    file << prop;
    file.close();
  }
}


void print_files(const File_spectrum_t& files, const HadSpec_t& spec, bool first)
{
  cout << __func__ << ": spec.had.size = " << spec.had.size() << endl;

  for(int i = 0; i < spec.had.size(); ++i){
    // mesons
    cout << __func__ << ": meson = " << files.had[i].meson.mes1.filename << endl;
    cout << __func__ << ": meson = " << files.had[i].meson.mes2.filename << endl;
    print_file(files.had[i].meson.mes1.filename, 
	       spec.had[i].meson.mes1, files.nprop, spec, first);
    print_file(files.had[i].meson.mes2.filename, 
	       spec.had[i].meson.mes2, files.nprop, spec, first);
    
    // baryons
    print_file(files.had[i].baryon.Lambda.filename,
	       spec.had[i].baryon.Lambda, files.nprop, spec, first);
    
    print_file(files.had[i].baryon.SigmaX.filename,
	       spec.had[i].baryon.SigmaX, files.nprop, spec, first);
    print_file(files.had[i].baryon.SigmaY.filename,
	       spec.had[i].baryon.SigmaY, files.nprop, spec, first);
    print_file(files.had[i].baryon.SigmaZ.filename,
	       spec.had[i].baryon.SigmaZ, files.nprop, spec, first);
  }
}


void construct_filenames(File_spectrum_t& files, const HadSpec_t& spec,
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

  for(int i = 0; i < spec.had.size(); ++i) 
  {
    const Wilson_hadron_measurements_t& had = spec.had[i];

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
      source_string = ".P";
      source_suffix = "P";
    }
    else if (had.Wl_src)
    {
      source_string = ".W";
      source_suffix = "W";
    }
    else if (had.Sl_src)
    {
      // Width names. Do not want trailing zeros
      bool sameP = true;
      Real wvf_param = had.forward_prop_headers.source_header_1.sourceSmearParam.wvf_param;
      {
	if (! ((had.forward_prop_headers.source_header_1.sourceSmearParam.wvf_param == wvf_param) &&
	       (had.forward_prop_headers.source_header_2.sourceSmearParam.wvf_param == wvf_param)))
	{
	  cout << "Not all source smearing widths are the same - cannot deal with that" << endl;
	  exit(1);
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

	source_string = ".D" + src_wvf_param_s;
	source_suffix = "S";
      }
      else
      {
	// This is not handled in a clean way - hack for now
	ostringstream wvf_param_str_1;
	wvf_param_str_1 << "G";
	wvf_param_str_1 << had.forward_prop_headers.source_header_1.sourceSmearParam.wvf_param;
	string wvf_ps_1 = wvf_param_str_1.str();

	string::size_type idx_1 = wvf_ps_1.find('.');  // replace '.' with 'p'
	if (idx_1 != string::npos)
	  wvf_ps_1[idx_1] = 'p';

	ostringstream wvf_param_str_2;
	wvf_param_str_2 << "G";
	wvf_param_str_2 << had.forward_prop_headers.source_header_1.sourceSmearParam.wvf_param;
	string wvf_ps_2 = wvf_param_str_2.str();

	string::size_type idx_2 = wvf_ps_2.find('.');  // replace '.' with 'p'
	if (idx_2 != string::npos)
	  wvf_ps_2[idx_2] = 'p';

	source_string = ".H" + wvf_ps_1 + "-" + wvf_ps_2;
	source_suffix = "S";
      }
    }
    else
    {
      cout << __func__ << ": Unsupported source type" << endl;
      exit(1);
    }

    cout << "source_string = " << source_string << endl;
      
    
    /*
     * Sink case
     */
    if (had.Pt_snk)
    {
      sink_string = ".P";
      sink_suffix = "P";
    }
    else if (had.Wl_snk)
    {
      sink_string = ".W";
      sink_suffix = "W";
    }
    else if (had.Sl_snk)
    {
      // Width names. Do not want trailing zeros
      ostringstream wvf_param_str;

      Real wvf_param = had.forward_prop_headers.sink_header_1.sinkSmearParam.wvf_param;
      {
	if (! ((had.forward_prop_headers.sink_header_1.sinkSmearParam.wvf_param == wvf_param) &&
	       (had.forward_prop_headers.sink_header_2.sinkSmearParam.wvf_param == wvf_param)))
	{
	  cout << "Not all sink smearing widths are the same - cannot deal with that" << endl;
	  cout << "snk_1 = " << had.forward_prop_headers.sink_header_1.sinkSmearParam.wvf_param << endl;
	  cout << "snk_2 = " << had.forward_prop_headers.sink_header_2.sinkSmearParam.wvf_param << endl;
	  exit(1);
	}
      }

      wvf_param_str << "G";
      wvf_param_str << wvf_param;
      string wvf_ps = wvf_param_str.str();

      string::size_type idx = wvf_ps.find('.');  // replace '.' with 'p'
      if (idx != string::npos)
	wvf_ps[idx] = 'p';

      snk_wvf_param_s = wvf_ps;

      sink_string = ".D" + snk_wvf_param_s;
      sink_suffix = "S";
    }
    else
    {
      cout << __func__ << ": Unsupported sink type" << endl;
      exit(1);
    }

    cout << "sink_string = " << sink_string << endl;

    cout << "  Mass_1 = " << Mass_1 << endl
	 << "  Mass_2 = " << Mass_2 << endl;

    cout << "source_string = " << source_string << endl;
    cout << "sink_string = " << sink_string << endl;
    
    //
    // If masses are degenerate, set a flag
    //
    bool degenerate_massP = (Mass_1 == Mass_2) ? true : false;


    //
    // Construct meson names
    //
    if (spec.param.MesonP)
    {

      /*
       * Create Meson source-sink smearing names
       */
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

	files.had[i].meson.mes1.filename = "hlmes1" 
	  + Mass_12 + source_string + sink_string  + "." + source_suffix + sink_suffix;
	files.had[i].meson.mes2.filename = "hlmes2" 
	  + Mass_12 + source_string + sink_string  + "." + source_suffix + sink_suffix;
	
	
      }
    }

    //
    // Construct baryon names
    //
    if (spec.param.BaryonP)
    {
      //if (degenerate_massP)
     

      /*
       * Create Baryon source-sink smearing names
       */
      cout << "Creating baryon file names" << endl;

      {
	string Mass_123;
	if (Mass_1 == Mass_2)
	{
	  Mass_123 = ".D" + Mass_1;
	}
	else
	{
	  Mass_123 = ".H" + Mass_1+"_"+Mass_2;
	}

	files.had[i].baryon.Lambda.filename = "LambdaQ"
	      + Mass_123 + source_string + sink_string + "." + source_suffix + sink_suffix;

	files.had[i].baryon.SigmaX.filename = "SigmaQx"
	  + Mass_123 + source_string + sink_string + "." + source_suffix + sink_suffix;
	files.had[i].baryon.SigmaY.filename = "SigmaQy"
	  + Mass_123 + source_string + sink_string + "." + source_suffix + sink_suffix;
	files.had[i].baryon.SigmaZ.filename = "SigmaQz"
	  + Mass_123 + source_string + sink_string + "." + source_suffix + sink_suffix;
      }
    }
  }

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

  string xml_in_root = "/StaticLightSpec";


  // Process the first file specially - read alls it params
  cout << "Open file " << argv[1] << endl;
  XMLReader xml_in(argv[1]);
  cout << "Done Open file " << argv[1] << endl;
  
  // Big nested structure that is image of entire file
  HadSpec_t  spec;

  // Read data
  cout << "Read config 1 data " << argv[1] << endl;
  read(xml_in, xml_in_root, spec);
  xml_in.close();
  
  cout << "spec.had.size = " << spec.had.size() << endl;
//  cout << "spec.had.meson.mom.size = " << spec.had[0].meson[0].mom.size() << endl;
//  cout << "spec.had.meson.mom.mesprop.size = " << spec.had[0].meson[0].mom[0].mesprop.size() << endl;
  

  // We care about the output version
  switch (spec.output_version.out_version)
  {
  case 14:
  case 15:
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

  exit(0);
}
