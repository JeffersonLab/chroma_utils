// $Id: strip_spectrum.cc,v 2.0 2008/12/05 04:44:04 edwards Exp $

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
  bool         MesonP;
  bool         CurrentP;
  bool         BaryonP;
  bool         HybMesP;
  bool         time_rev;
  int          mom2_max;
  bool         avg_equiv_mom;

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


struct Meson_momenta_t
{
  int sink_mom_num;
  Array<int> sink_mom;
  Array<Real> mesprop;
};

struct Wilson_hybmes_t
{
  int kv;
  Array<Meson_momenta_t> mom;
};

struct Wilson_mesons_t
{
  int gamma_value;
  Array<Meson_momenta_t> mom;
};

struct Vector_current_t
{
  Array<Real>    vector;
};

struct Axial_current_t
{
  Array<Real>    axial;
};

struct Wilson_currents_t
{
  int                      num_vec_cur;
  Array<Vector_current_t>  vector;
  Array<Axial_current_t>   axial;
};

struct Baryon_momenta_t
{
  int            sink_mom_num;
  Array<int>     sink_mom;
  Array<Complex> barprop;
};

struct Wilson_baryons_t
{
  int                     baryon_num;
  Array<Baryon_momenta_t> mom;
};

struct Wilson_hadron_measurements_t
{
  int loop;

  PropSource_t            source_header;
  Real                    Mass;       // quark mass (bare units)

  bool         Pt_src;   // Not read, but determined from source
  bool         Sl_src;   // Not read, but determined from source
  bool         Wl_src;   // Not read, but determined from source

  Array<Wilson_mesons_t>  meson_PP;
  Array<Wilson_mesons_t>  meson_PS;
  Array<Wilson_mesons_t>  meson_PW;
  Array<Wilson_mesons_t>  meson_SP;
  Array<Wilson_mesons_t>  meson_SS;
  Array<Wilson_mesons_t>  meson_SW;
  Array<Wilson_mesons_t>  meson_WP;
  Array<Wilson_mesons_t>  meson_WS;
  Array<Wilson_mesons_t>  meson_WW;

  Array<Wilson_hybmes_t>  hybmes_PP;
  Array<Wilson_hybmes_t>  hybmes_PS;
  Array<Wilson_hybmes_t>  hybmes_PW;
  Array<Wilson_hybmes_t>  hybmes_SP;
  Array<Wilson_hybmes_t>  hybmes_SS;
  Array<Wilson_hybmes_t>  hybmes_SW;
  Array<Wilson_hybmes_t>  hybmes_WP;
  Array<Wilson_hybmes_t>  hybmes_WS;
  Array<Wilson_hybmes_t>  hybmes_WW;

  Wilson_currents_t       current_PP;
  Wilson_currents_t       current_SP;
  Wilson_currents_t       current_WP;

  Array<Wilson_baryons_t> baryon_PP;
  Array<Wilson_baryons_t> baryon_PS;
  Array<Wilson_baryons_t> baryon_PW;
  Array<Wilson_baryons_t> baryon_SP;
  Array<Wilson_baryons_t> baryon_SS;
  Array<Wilson_baryons_t> baryon_SW;
  Array<Wilson_baryons_t> baryon_WP;
  Array<Wilson_baryons_t> baryon_WS;
  Array<Wilson_baryons_t> baryon_WW;

};

struct Spectrum_w_t
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
  param.HybMesP = false;

  switch (version)
  {
  case 9:
    param.Wl_snk = false;
    break;

  case 10:
    read(top, "Wl_snk", param.Wl_snk);
    break;

  case 11:
  case 12:
  case 13:
    // At the moment, the new hybrids stuff is not supported.
    read(top, "Wl_snk", param.Wl_snk);
    read(top, "HybMesP", param.HybMesP);
    break;

  default:
    error_stream << "Within Param read - unsupported param version = "
		 << version << endl;
    throw error_stream.str();
  }

  read(top, "Pt_snk", param.Pt_snk);
  read(top, "Sl_snk", param.Sl_snk);
  read(top, "MesonP", param.MesonP);
  read(top, "CurrentP", param.CurrentP);
  read(top, "BaryonP", param.BaryonP);
  read(top, "time_rev", param.time_rev);
  read(top, "wvf_kind", param.wvf_kind);
  read(top, "wvf_param", param.wvf_param);
  read(top, "wvfIntPar", param.wvfIntPar);

  if (param.wvf_param.size() != param.wvfIntPar.size())
  {
    error_stream << "Unexpected wvf_param array size";
    throw error_stream.str();
  }

  read(top, "mom2_max", param.mom2_max);
  read(top, "avg_equiv_mom", param.avg_equiv_mom);

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

  read(top, "gamma_value", mes.gamma_value);
  read(top, "momenta", mes.mom);
}

// hybrid meson struct
void read(XMLReader& xml, const string& path, Wilson_hybmes_t& mes)
{
  XMLReader top(xml, path);

  read(top, "kv", mes.kv);
  read(top, "momenta", mes.mom);
}

// current struct
void read(XMLReader& xml, const string& path, Vector_current_t& cur)
{
  XMLReader top(xml, path);

  read(top, "vector_current", cur.vector);
}

// current struct
void read(XMLReader& xml, const string& path, Axial_current_t& cur)
{
  XMLReader top(xml, path);

  read(top, "axial_current", cur.axial);
}

// current struct
void read(XMLReader& xml, const string& path, Wilson_currents_t& cur)
{
  XMLReader top(xml, path);

  string xpath;

  // This should be here in version 7 and above
  xpath = "Vector_currents";
  if (top.count(xpath) != 0)
    read(top, xpath, cur.vector);

  xpath = "Axial_currents";
  if (top.count(xpath) != 0)
    read(top, xpath, cur.axial);
}

// baryon momenta
void read(XMLReader& xml, const string& path, Baryon_momenta_t& mom)
{
  XMLReader top(xml, path);

  read(top, "sink_mom_num", mom.sink_mom_num);
  read(top, "sink_mom", mom.sink_mom);
  read(top, "barprop", mom.barprop);
}

// baryon struct
void read(XMLReader& xml, const string& path, Wilson_baryons_t& bar)
{
  XMLReader top(xml, path);

  read(top, "baryon_num", bar.baryon_num);
  read(top, "momenta", bar.mom);
}

// Read a hadron measurements
void read(XMLReader& xml, const string& path, Wilson_hadron_measurements_t& had)
{
//  cout << "Wilson_had = " << path << endl;

  XMLReader top(xml, path);

  read(top, "loop", had.loop);

  read(top, "PropSource", had.source_header);
  if (top.count("Mass_mes") != 0)
  {
    read(top, "Mass_mes", had.Mass);
  }
  else
  {
    ChromaProp_t prop_header;
    read(top, "ForwardProp", prop_header);
    had.Mass = prop_header.Mass;
  }

  // Determine what kind of source to use
  had.Pt_src = (had.source_header.source_type == "POINT_SOURCE") ? true : false;
  had.Sl_src = (had.source_header.source_type == "SHELL_SOURCE") ? true : false;
  had.Wl_src = (had.source_header.source_type == "WALL_SOURCE")  ? true : false;

  string xpath;

  // Mesons
  xpath = "Point_Point_Wilson_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.meson_PP);

  xpath = "Point_Shell_Wilson_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.meson_PS);

  xpath = "Point_Wall_Wilson_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.meson_PW);

  xpath = "Shell_Point_Wilson_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.meson_SP);

  xpath = "Shell_Shell_Wilson_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.meson_SS);

  xpath = "Shell_Wall_Wilson_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.meson_SW);

  xpath = "Wall_Point_Wilson_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.meson_WP);

  xpath = "Wall_Wall_Wilson_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.meson_WS);

  xpath = "Wall_Wall_Wilson_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.meson_WW);

  // Currents
  xpath = "Point_Point_Meson_Currents";
  if (top.count(xpath) != 0)
    read(top, xpath, had.current_PP);

  xpath = "Shell_Point_Meson_Currents";
  if (top.count(xpath) != 0)
    read(top, xpath, had.current_SP);

  xpath = "Wall_Point_Meson_Currents";
  if (top.count(xpath) != 0)
    read(top, xpath, had.current_WP);

  // Baryons
  xpath = "Point_Point_Wilson_Baryons"; 
  if (top.count(xpath) != 0)
    read(top, xpath, had.baryon_PP);

  xpath = "Point_Shell_Wilson_Baryons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.baryon_PS);

  xpath = "Point_Wall_Wilson_Baryons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.baryon_PW);

  xpath = "Shell_Point_Wilson_Baryons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.baryon_SP);

  xpath = "Shell_Shell_Wilson_Baryons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.baryon_SS);

  xpath = "Shell_Wall_Wilson_Baryons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.baryon_SW);

  xpath = "Wall_Point_Wilson_Baryons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.baryon_WP);

  xpath = "Wall_Shell_Wilson_Baryons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.baryon_WS);

  xpath = "Wall_Wall_Wilson_Baryons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.baryon_WW);

  // Hybrid mesons
  xpath = "Point_Point_Wilson_Hybrid_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.hybmes_PP);

  xpath = "Point_Shell_Wilson_Hybrid_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.hybmes_PS);

  xpath = "Point_Wall_Wilson_Hybrid_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.hybmes_PW);

  xpath = "Shell_Point_Wilson_Hybrid_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.hybmes_SP);

  xpath = "Shell_Shell_Wilson_Hybrid_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.hybmes_SS);

  xpath = "Shell_Wall_Wilson_Hybrid_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.hybmes_SW);

  xpath = "Wall_Point_Wilson_Hybrid_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.hybmes_WP);

  xpath = "Wall_Wall_Wilson_Hybrid_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.hybmes_WS);

  xpath = "Wall_Wall_Wilson_Hybrid_Mesons";
  if (top.count(xpath) != 0)
    read(top, xpath, had.hybmes_WW);
}

// Mother of all readers
void read(XMLReader& xml, const string& path, Spectrum_w_t& spec)
{
  ostringstream error_stream;

  try 
  {
    XMLReader mother(xml, path);

    read(mother, "Output_version", spec.output_version);

    read(mother, "Output_version", spec.output_version);
    if (mother.count("Input/spectrum_w/Param") != 0)
      read(mother, "Input/spectrum_w/Param", spec.param);
    else if (mother.count("Input/Param") != 0)
      read(mother, "Input/Param", spec.param);
    else
    {
      throw string("Cannot find Param input");
    }

    switch (spec.output_version.out_version)
    {
    case 9:
    case 10:
    case 11:
    case 12:
    case 13:
      break;

    default:
      error_stream << "read(spectrum): Unsupported output version " 
		   << spec.output_version.out_version << endl;
      throw error_stream.str();
    }

    read(mother, "Wilson_hadron_measurements", spec.had);
  }
  catch( const string& error) { 
    cerr << "Error reading Spectrum_w_t: " << error << endl;
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
  Array<File_t> mom;
};

struct File_current_t
{
  Array<File_t> vector;
  Array<File_t> axial;
};

struct File_hadron_measurements_t
{
  Array<File_particles_t> meson_PP;
  Array<File_particles_t> meson_PS;
  Array<File_particles_t> meson_PW;
  Array<File_particles_t> meson_SP;
  Array<File_particles_t> meson_SS;
  Array<File_particles_t> meson_SW;
  Array<File_particles_t> meson_WP;
  Array<File_particles_t> meson_WS;
  Array<File_particles_t> meson_WW;

  Array<File_particles_t> hybmes_PP;
  Array<File_particles_t> hybmes_PS;
  Array<File_particles_t> hybmes_PW;
  Array<File_particles_t> hybmes_SP;
  Array<File_particles_t> hybmes_SS;
  Array<File_particles_t> hybmes_SW;
  Array<File_particles_t> hybmes_WP;
  Array<File_particles_t> hybmes_WS;
  Array<File_particles_t> hybmes_WW;

  File_current_t          current_PP;
  File_current_t          current_SP;
  File_current_t          current_WP;

  Array<File_particles_t> baryon_PP;
  Array<File_particles_t> baryon_PS;
  Array<File_particles_t> baryon_PW;
  Array<File_particles_t> baryon_SP;
  Array<File_particles_t> baryon_SS;
  Array<File_particles_t> baryon_SW;
  Array<File_particles_t> baryon_WP;
  Array<File_particles_t> baryon_WS;
  Array<File_particles_t> baryon_WW;
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


string mom_name(const Array<int>& sink_mom)
{
  ostringstream s;

  if (! (sink_mom[0] == 0 && sink_mom[1] == 0 && sink_mom[2] == 0) )
    s << "_px" << sink_mom[0] << "_py" << sink_mom[1] << "_pz" << sink_mom[2];

  return s.str();
}


void print_file(const string& filename, const Array<Real>& mesprop, 
		int nprop, 
		const Spectrum_w_t& spec, 
		bool first)
{
  int j_decay = spec.had[0].source_header.j_decay;

  if (first)
  {
    ofstream file(filename.c_str());

    // Write header line
    file << nprop << " " << spec.param.nrow[j_decay] << " 0 " 
	 << spec.param.nrow[0] << " 1" << endl;

    file << mesprop;   // Write the file
    file.close();
  }
  else
  {
    ofstream file(filename.c_str(), std::ios_base::app);
    file << mesprop;
    file.close();
  }
}


// Complex writer for now only writes real part   !!!!
void print_file(const string& filename, const Array<Complex>& barprop, 
		int nprop, 
		const Spectrum_w_t& spec, 
		bool first)
{
  int j_decay = spec.had[0].source_header.j_decay;

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
      // mesons
      for(int m = 0; m < spec.had[i].meson_PP.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].meson_PP[m].mom.size(); ++n)
	  print_file(files.had[i].meson_PP[m].mom[n].filename, 
		     spec.had[i].meson_PP[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // hybrid mesons
      for(int m = 0; m < spec.had[i].hybmes_PP.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].hybmes_PP[m].mom.size(); ++n)
	  print_file(files.had[i].hybmes_PP[m].mom[n].filename, 
		     spec.had[i].hybmes_PP[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // vector currents
      for(int n = 0; n < spec.had[i].current_PP.vector.size(); ++n)
      {
	print_file(files.had[i].current_PP.vector[n].filename, 
		   spec.had[i].current_PP.vector[n].vector, 
		   files.nprop, spec, first);
      }

      // axial currents
      for(int n = 0; n < spec.had[i].current_PP.axial.size(); ++n)
      {
	print_file(files.had[i].current_PP.axial[n].filename, 
		   spec.had[i].current_PP.axial[n].axial, 
		   files.nprop, spec, first);
      }

      // baryons
      for(int m = 0; m < spec.had[i].baryon_PP.size(); ++m)
	for(int n = 0; n < spec.had[i].baryon_PP[m].mom.size(); ++n)
	  print_file(files.had[i].baryon_PP[m].mom[n].filename,
		     spec.had[i].baryon_PP[m].mom[n].barprop, 
		     files.nprop, spec, first);
    }

    if ( spec.had[i].Pt_src && spec.param.Sl_snk ) 
    {
      // mesons
      for(int m = 0; m < spec.had[i].meson_PS.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].meson_PS[m].mom.size(); ++n)
	  print_file(files.had[i].meson_PS[m].mom[n].filename,
		     spec.had[i].meson_PS[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // hybrid mesons
      for(int m = 0; m < spec.had[i].hybmes_PS.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].hybmes_PS[m].mom.size(); ++n)
	  print_file(files.had[i].hybmes_PS[m].mom[n].filename,
		     spec.had[i].hybmes_PS[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // baryons
      for(int m = 0; m < spec.had[i].baryon_PS.size(); ++m)
	for(int n = 0; n < spec.had[i].baryon_PS[m].mom.size(); ++n)
	  print_file(files.had[i].baryon_PS[m].mom[n].filename,
		     spec.had[i].baryon_PS[m].mom[n].barprop, 
		     files.nprop, spec, first);
    }

    if ( spec.had[i].Pt_src && spec.param.Wl_snk ) 
    {
      // mesons
      for(int m = 0; m < spec.had[i].meson_PW.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].meson_PW[m].mom.size(); ++n)
	  print_file(files.had[i].meson_PW[m].mom[n].filename,
		     spec.had[i].meson_PW[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // hybrid mesons
      for(int m = 0; m < spec.had[i].hybmes_PW.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].hybmes_PW[m].mom.size(); ++n)
	  print_file(files.had[i].hybmes_PW[m].mom[n].filename,
		     spec.had[i].hybmes_PW[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // baryons
      for(int m = 0; m < spec.had[i].baryon_PW.size(); ++m)
	for(int n = 0; n < spec.had[i].baryon_PW[m].mom.size(); ++n)
	  print_file(files.had[i].baryon_PW[m].mom[n].filename,
		     spec.had[i].baryon_PW[m].mom[n].barprop, 
		     files.nprop, spec, first);
    }

    if ( spec.had[i].Sl_src && spec.param.Pt_snk ) 
    {
      // mesons
      for(int m = 0; m < spec.had[i].meson_SP.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].meson_SP[m].mom.size(); ++n)
	  print_file(files.had[i].meson_SP[m].mom[n].filename,
		     spec.had[i].meson_SP[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // hybrid mesons
      for(int m = 0; m < spec.had[i].hybmes_SP.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].hybmes_SP[m].mom.size(); ++n)
	  print_file(files.had[i].hybmes_SP[m].mom[n].filename,
		     spec.had[i].hybmes_SP[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // vector currents
      for(int n = 0; n < spec.had[i].current_SP.vector.size(); ++n)
      {
	print_file(files.had[i].current_SP.vector[n].filename, 
		   spec.had[i].current_SP.vector[n].vector, 
		   files.nprop, spec, first);
      }

      // axial currents
      for(int n = 0; n < spec.had[i].current_SP.axial.size(); ++n)
      {
	print_file(files.had[i].current_SP.axial[n].filename, 
		   spec.had[i].current_SP.axial[n].axial, 
		   files.nprop, spec, first);
      }

      // baryons
      for(int m = 0; m < spec.had[i].baryon_SP.size(); ++m)
	for(int n = 0; n < spec.had[i].baryon_SP[m].mom.size(); ++n)
	  print_file(files.had[i].baryon_SP[m].mom[n].filename,
		     spec.had[i].baryon_SP[m].mom[n].barprop, 
		     files.nprop, spec, first);
    }
    
    if ( spec.had[i].Sl_src && spec.param.Sl_snk ) 
    {
      // mesons
      for(int m = 0; m < spec.had[i].meson_SS.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].meson_SS[m].mom.size(); ++n)
	  print_file(files.had[i].meson_SS[m].mom[n].filename,
		     spec.had[i].meson_SS[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // hybrid mesons
      for(int m = 0; m < spec.had[i].hybmes_SS.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].hybmes_SS[m].mom.size(); ++n)
	  print_file(files.had[i].hybmes_SS[m].mom[n].filename,
		     spec.had[i].hybmes_SS[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // baryons
      for(int m = 0; m < spec.had[i].baryon_SS.size(); ++m)
	for(int n = 0; n < spec.had[i].baryon_SS[m].mom.size(); ++n)
	  print_file(files.had[i].baryon_SS[m].mom[n].filename,
		     spec.had[i].baryon_SS[m].mom[n].barprop, 
		     files.nprop, spec, first);
    } 
    
    if ( spec.had[i].Sl_src && spec.param.Wl_snk ) 
    {
      // mesons
      for(int m = 0; m < spec.had[i].meson_SW.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].meson_SW[m].mom.size(); ++n)
	  print_file(files.had[i].meson_SW[m].mom[n].filename,
		     spec.had[i].meson_SW[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // hybrid mesons
      for(int m = 0; m < spec.had[i].hybmes_SW.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].hybmes_SW[m].mom.size(); ++n)
	  print_file(files.had[i].hybmes_SW[m].mom[n].filename,
		     spec.had[i].hybmes_SW[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // baryons
      for(int m = 0; m < spec.had[i].baryon_SW.size(); ++m)
	for(int n = 0; n < spec.had[i].baryon_SW[m].mom.size(); ++n)
	  print_file(files.had[i].baryon_SW[m].mom[n].filename,
		     spec.had[i].baryon_SW[m].mom[n].barprop, 
		     files.nprop, spec, first);
    } 

    if ( spec.had[i].Wl_src && spec.param.Pt_snk ) 
    {
      // mesons
      for(int m = 0; m < spec.had[i].meson_WP.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].meson_WP[m].mom.size(); ++n)
	  print_file(files.had[i].meson_WP[m].mom[n].filename,
		     spec.had[i].meson_WP[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // hybrid mesons
      for(int m = 0; m < spec.had[i].hybmes_WP.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].hybmes_WP[m].mom.size(); ++n)
	  print_file(files.had[i].hybmes_WP[m].mom[n].filename,
		     spec.had[i].hybmes_WP[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // vector currents
      for(int n = 0; n < spec.had[i].current_WP.vector.size(); ++n)
      {
	print_file(files.had[i].current_WP.vector[n].filename, 
		   spec.had[i].current_WP.vector[n].vector, 
		   files.nprop, spec, first);
      }

      // axial currents
      for(int n = 0; n < spec.had[i].current_WP.axial.size(); ++n)
      {
	print_file(files.had[i].current_WP.axial[n].filename, 
		   spec.had[i].current_WP.axial[n].axial, 
		   files.nprop, spec, first);
      }

      // baryons
      for(int m = 0; m < spec.had[i].baryon_WP.size(); ++m)
	for(int n = 0; n < spec.had[i].baryon_WP[m].mom.size(); ++n)
	  print_file(files.had[i].baryon_WP[m].mom[n].filename,
		     spec.had[i].baryon_WP[m].mom[n].barprop, 
		     files.nprop, spec, first);
    }
    
    if ( spec.had[i].Wl_src && spec.param.Sl_snk ) 
    {
      // mesons
      for(int m = 0; m < spec.had[i].meson_WS.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].meson_WS[m].mom.size(); ++n)
	  print_file(files.had[i].meson_WS[m].mom[n].filename,
		     spec.had[i].meson_WS[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // hybrid mesons
      for(int m = 0; m < spec.had[i].hybmes_WS.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].hybmes_WS[m].mom.size(); ++n)
	  print_file(files.had[i].hybmes_WS[m].mom[n].filename,
		     spec.had[i].hybmes_WS[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // baryons
      for(int m = 0; m < spec.had[i].baryon_WS.size(); ++m)
	for(int n = 0; n < spec.had[i].baryon_WS[m].mom.size(); ++n)
	  print_file(files.had[i].baryon_WS[m].mom[n].filename,
		     spec.had[i].baryon_WS[m].mom[n].barprop, 
		     files.nprop, spec, first);
    } 
    
    if ( spec.had[i].Wl_src && spec.param.Wl_snk ) 
    {
      // mesons
      for(int m = 0; m < spec.had[i].meson_WW.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].meson_WW[m].mom.size(); ++n)
	  print_file(files.had[i].meson_WW[m].mom[n].filename,
		     spec.had[i].meson_WW[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // hybrid mesons
      for(int m = 0; m < spec.had[i].hybmes_WW.size(); ++m)
      {
	for(int n = 0; n < spec.had[i].hybmes_WW[m].mom.size(); ++n)
	  print_file(files.had[i].hybmes_WW[m].mom[n].filename,
		     spec.had[i].hybmes_WW[m].mom[n].mesprop, 
		     files.nprop, spec, first);
      }

      // baryons
      for(int m = 0; m < spec.had[i].baryon_WW.size(); ++m)
	for(int n = 0; n < spec.had[i].baryon_WW[m].mom.size(); ++n)
	  print_file(files.had[i].baryon_WW[m].mom[n].filename,
		     spec.had[i].baryon_WW[m].mom[n].barprop, 
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
    Mass_str << int(10000*spec.had[i].Mass + 0.5);
    Mass_s[i] = Mass_str.str();

    ostringstream error_stream;

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
	error_stream << "Unknown or unsupported source wvf_kind";
	throw error_stream.str();
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
	error_stream << "Unknown or unsupported sink wvf_kind" ;
	throw error_stream.str();
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


  //
  // Construct meson names
  //
  if (spec.param.MesonP)
  {
    Array<string> meson_particle(Ns*Ns);    // Ns*Ns

    meson_particle[0]  = "a0";
    meson_particle[1]  = "rho_x";
    meson_particle[2]  = "rho_y";
    meson_particle[3]  = "b1_z";
    meson_particle[4]  = "rho_z";
    meson_particle[5]  = "b1_y";
    meson_particle[6]  = "b1_x";
    meson_particle[7]  = "pion";
    meson_particle[8]  = "a0";
    meson_particle[9]  = "rho_x";
    meson_particle[10] = "rho_y";
    meson_particle[11] = "a1_z";
    meson_particle[12] = "rho_z";
    meson_particle[13] = "a1_y";
    meson_particle[14] = "a1_x";
    meson_particle[15] = "pion";

    Array<string> meson_smear_state(16);    // Ns*Ns

    meson_smear_state[0]  = "_1" ;
    meson_smear_state[1]  = "_1" ;
    meson_smear_state[2]  = "_1" ;
    meson_smear_state[3]  = "_1" ;
    meson_smear_state[4]  = "_1" ;
    meson_smear_state[5]  = "_1" ;
    meson_smear_state[6]  = "_1" ;
    meson_smear_state[7]  = "_2" ;
    meson_smear_state[8]  = "_2" ;
    meson_smear_state[9]  = "_2" ;
    meson_smear_state[10] = "_2" ;
    meson_smear_state[11] = "_1" ;
    meson_smear_state[12] = "_2" ;
    meson_smear_state[13] = "_1" ;
    meson_smear_state[14] = "_1" ;
    meson_smear_state[15] = "_1" ;

    /*
     * Create Meson source-sink smearing names
     */
    cout << "Creating meson file names" << endl;

    // For now assume only diagonal valence mesons
    for(int i = 0; i < spec.had.size(); ++i) 
    {
      int i1=i, i2=i;
      bool diag = true;
      string Mass_12 = ".D" + Mass_s[i1];
      string wvf_param12_s = "";   // ignore for the moment, should fill in

      if ( spec.had[i].Pt_src && spec.param.Pt_snk ) 
      {
	files.had[i].meson_PP.resize(spec.had[i].meson_PP.size());

	for(int m=0; m < spec.had[i].meson_PP.size(); ++m)
	{
	  string ms = meson_smear_state[m];
	
	  files.had[i].meson_PP[m].mom.resize(spec.had[i].meson_PP[m].mom.size());
	
	  for(int n = 0; n < spec.had[i].meson_PP[m].mom.size(); ++n)
	  {
	    files.had[i].meson_PP[m].mom[n].filename = meson_particle[m] 
	      + mom_name(spec.had[i].meson_PP[m].mom[n].sink_mom) 
	      + Mass_12 + ".P" + ms + ".P" + ms + ".PP";
	  }
	}
      }

      if ( diag )
      {
	if ( spec.had[i].Pt_src && spec.param.Sl_snk ) 
	{
	  files.had[i].meson_PS.resize(spec.had[i].meson_PS.size());

	  for(int m=0; m < spec.had[i].meson_PS.size(); ++m)
	  {
	    string ms = meson_smear_state[m];

	    files.had[i].meson_PS[m].mom.resize(spec.had[i].meson_PS[m].mom.size());

	    for(int n = 0; n < spec.had[i].meson_PS[m].mom.size(); ++n)
	    {
	      files.had[i].meson_PS[m].mom[n].filename = meson_particle[m] 
		+ mom_name(spec.had[i].meson_PS[m].mom[n].sink_mom) 
		+ Mass_12 + ".P" + ms + ".D" + snk_wvf_param_s[i2] + ms + ".PS";
	    }
	  }
	}

	if ( spec.had[i].Pt_src && spec.param.Wl_snk ) 
	{
	  files.had[i].meson_PW.resize(spec.had[i].meson_PW.size());

	  for(int m=0; m < spec.had[i].meson_PW.size(); ++m)
	  {
	    string ms = meson_smear_state[m];
	
	    files.had[i].meson_PW[m].mom.resize(spec.had[i].meson_PW[m].mom.size());
	
	    for(int n = 0; n < spec.had[i].meson_PW[m].mom.size(); ++n)
	    {
	      files.had[i].meson_PW[m].mom[n].filename = meson_particle[m] 
		+ mom_name(spec.had[i].meson_PW[m].mom[n].sink_mom) 
		+ Mass_12 + ".P" + ms + ".W" + ms + ".PW";
	    }
	  }
	}

	if ( spec.had[i].Sl_src && spec.param.Pt_snk ) 
	{
	  files.had[i].meson_SP.resize(spec.had[i].meson_SP.size());

	  for(int m=0; m < spec.had[i].meson_SP.size(); ++m)
	  {
	    string ms = meson_smear_state[m];

	    files.had[i].meson_SP[m].mom.resize(spec.had[i].meson_SP[m].mom.size());

	    for(int n = 0; n < spec.had[i].meson_SP[m].mom.size(); ++n)
	    {
	      files.had[i].meson_SP[m].mom[n].filename = meson_particle[m] 
		+ mom_name(spec.had[i].meson_SP[m].mom[n].sink_mom) 
		+ Mass_12 + ".D" + src_wvf_param_s[i2] + ms + ".P" + ms + ".SP";

//	    cout << spec.had[i].meson_SP[m].mom[n].filename << endl;
	    }
	  }
	}

	if ( spec.had[i].Sl_src && spec.param.Sl_snk ) 
	{
	  files.had[i].meson_SS.resize(spec.had[i].meson_SS.size());

	  for(int m=0; m < spec.had[i].meson_SS.size(); ++m)
	  {
	    string ms = meson_smear_state[m];

	    files.had[i].meson_SS[m].mom.resize(spec.had[i].meson_SS[m].mom.size());

	    for(int n = 0; n < spec.had[i].meson_SS[m].mom.size(); ++n)
	    {
	      files.had[i].meson_SS[m].mom[n].filename = meson_particle[m] 
		+ mom_name(spec.had[i].meson_SS[m].mom[n].sink_mom) 
		+ Mass_12 + ".D" + src_wvf_param_s[i1] + ms + ".D" + snk_wvf_param_s[i2] + ms + ".SS";

//	    cout << spec.had[i].meson_SP[m].mom[n].filename << endl;
	    }
	  }
	}

	if ( spec.had[i].Sl_src && spec.param.Wl_snk ) 
	{
	  files.had[i].meson_SW.resize(spec.had[i].meson_SW.size());

	  for(int m=0; m < spec.had[i].meson_SW.size(); ++m)
	  {
	    string ms = meson_smear_state[m];

	    files.had[i].meson_SW[m].mom.resize(spec.had[i].meson_SW[m].mom.size());

	    for(int n = 0; n < spec.had[i].meson_SW[m].mom.size(); ++n)
	    {
	      files.had[i].meson_SW[m].mom[n].filename = meson_particle[m] 
		+ mom_name(spec.had[i].meson_SW[m].mom[n].sink_mom) 
		+ Mass_12 + ".D" + src_wvf_param_s[i2] + ms + ".W" + ms + ".SW";

//	    cout << spec.had[i].meson_SW[m].mom[n].filename << endl;
	    }
	  }
	}

	if ( spec.had[i].Wl_src && spec.param.Pt_snk ) 
	{
	  files.had[i].meson_WP.resize(spec.had[i].meson_WP.size());

	  for(int m=0; m < spec.had[i].meson_WP.size(); ++m)
	  {
	    string ms = meson_smear_state[m];
	
	    files.had[i].meson_WP[m].mom.resize(spec.had[i].meson_WP[m].mom.size());
	
	    for(int n = 0; n < spec.had[i].meson_WP[m].mom.size(); ++n)
	    {
	      files.had[i].meson_WP[m].mom[n].filename = meson_particle[m] 
		+ mom_name(spec.had[i].meson_WP[m].mom[n].sink_mom) 
		+ Mass_12 + ".W" + ms + ".P" + ms + ".WP";
	    }
	  }
	}

	if ( spec.had[i].Wl_src && spec.param.Sl_snk ) 
	{
	  files.had[i].meson_WS.resize(spec.had[i].meson_WS.size());

	  for(int m=0; m < spec.had[i].meson_WS.size(); ++m)
	  {
	    string ms = meson_smear_state[m];

	    files.had[i].meson_WS[m].mom.resize(spec.had[i].meson_WS[m].mom.size());

	    for(int n = 0; n < spec.had[i].meson_WS[m].mom.size(); ++n)
	    {
	      files.had[i].meson_WS[m].mom[n].filename = meson_particle[m] 
		+ mom_name(spec.had[i].meson_WS[m].mom[n].sink_mom) 
		+ Mass_12 + ".W" + ms + ".D" + snk_wvf_param_s[i2] + ms + ".WS";
	    }
	  }
	}

	if ( spec.had[i].Wl_src && spec.param.Wl_snk ) 
	{
	  files.had[i].meson_WW.resize(spec.had[i].meson_WW.size());

	  for(int m=0; m < spec.had[i].meson_WW.size(); ++m)
	  {
	    string ms = meson_smear_state[m];
	
	    files.had[i].meson_WW[m].mom.resize(spec.had[i].meson_WW[m].mom.size());
	
	    for(int n = 0; n < spec.had[i].meson_WW[m].mom.size(); ++n)
	    {
	      files.had[i].meson_WW[m].mom[n].filename = meson_particle[m] 
		+ mom_name(spec.had[i].meson_WW[m].mom[n].sink_mom) 
		+ Mass_12 + ".W" + ms + ".W" + ms + ".WW";
	    }
	  }
	}
      } 
      else    // if (! diag )   NOTE: the code will never get to this part for now
      {
	if ( spec.had[i].Pt_src && spec.param.Sl_snk ) 
	{
	  files.had[i].meson_PS.resize(spec.had[i].meson_PS.size());

	  for(int m=0; m < spec.had[i].meson_PS.size(); ++m)
	  {
	    string ms = meson_smear_state[m];

	    files.had[i].meson_PS[m].mom.resize(spec.had[i].meson_PS[m].mom.size());

	    for(int n = 0; n < spec.had[i].meson_PS[m].mom.size(); ++n)
	    {
	      files.had[i].meson_PS[m].mom[n].filename = meson_particle[m] 
		+ mom_name(spec.had[i].meson_PS[m].mom[n].sink_mom) 
		+ Mass_12 + ".P" + ms + wvf_param12_s + ms + ".PS";
	    }
	  }
	}

	if ( spec.had[i].Sl_src && spec.param.Pt_snk ) 
	{
	  files.had[i].meson_SP.resize(spec.had[i].meson_SP.size());

	  for(int m=0; m < spec.had[i].meson_SP.size(); ++m)
	  {
	    string ms = meson_smear_state[m];

	    files.had[i].meson_SP[m].mom.resize(spec.had[i].meson_SP[m].mom.size());

	    for(int n = 0; n < spec.had[i].meson_SP[m].mom.size(); ++n)
	    {
	      files.had[i].meson_SP[m].mom[n].filename = meson_particle[m] 
		+ mom_name(spec.had[i].meson_SP[m].mom[n].sink_mom) 
		+ Mass_12 + wvf_param12_s + ms + ".P" + ms + ".SP";
	    }
	  }
	}

	if ( spec.had[i].Sl_src && spec.param.Sl_snk ) 
	{
	  files.had[i].meson_SS.resize(spec.had[i].meson_SS.size());

	  for(int m=0; m < spec.had[i].meson_SS.size(); ++m)
	  {
	    string ms = meson_smear_state[m];

	    files.had[i].meson_SS[m].mom.resize(spec.had[i].meson_SS[m].mom.size());

	    for(int n = 0; n < spec.had[i].meson_SS[m].mom.size(); ++n)
	    {
	      files.had[i].meson_SS[m].mom[n].filename = meson_particle[m] 
		+ mom_name(spec.had[i].meson_SS[m].mom[n].sink_mom) 
		+ Mass_12 + wvf_param12_s + ms + wvf_param12_s + ms + ".SS";
	    }
	  }
	}

      }  // end if (diag == 1)
    }  // end for i
  }

 
  //
  // Construct current names
  //
  if (spec.param.CurrentP)
  {
    Array<string> vector_particle(12);

    vector_particle[0]  = "nonconserved_vector_x";
    vector_particle[1]  = "nonconserved_vector_y";
    vector_particle[2]  = "nonconserved_vector_z";
    vector_particle[3]  = "conserved_vector_x";
    vector_particle[4]  = "conserved_vector_y";
    vector_particle[5]  = "conserved_vector_z";
    vector_particle[6]  = "improved_rho_vector_x";
    vector_particle[7]  = "improved_rho_vector_y";
    vector_particle[8]  = "improved_rho_vector_z";
    vector_particle[9]  = "improved_rho_vector_x";
    vector_particle[10] = "improved_rho_vector_y";
    vector_particle[11] = "improved_rho_vector_z";

    Array<string> axial_particle(2);
    axial_particle[0] = "nonlocal_axial";
    axial_particle[1] = "local_axial";

    /*
     * Create Current source-sink smearing names
     */
    cout << "Creating current file names" << endl;

    // For now assume only diagonal valence currents
    for(int i = 0; i < spec.had.size(); ++i) 
    {
      int i1=i, i2=i;
      bool diag = true;
      string Mass_12 = ".D" + Mass_s[i1];
      string wvf_param12_s = "";   // ignore for the moment, should fill in

      if ( spec.had[i].Pt_src && spec.param.Pt_snk ) 
      {
	string source_sink = Mass_12 + ".P.P.PP";
	int nn;

	files.had[i].current_PP.vector.resize(spec.had[i].current_PP.vector.size());
	nn = spec.had[i].current_PP.vector.size();
	for(int n = 0; n < nn; ++n)
	  files.had[i].current_PP.vector[n].filename = vector_particle[n] + source_sink;

	files.had[i].current_PP.axial.resize(spec.had[i].current_PP.axial.size());
	nn = spec.had[i].current_PP.axial.size();
	for(int n = 0; n < nn; ++n)
	  files.had[i].current_PP.axial[n].filename = axial_particle[n] + source_sink;
      }

      if ( spec.had[i].Sl_src && spec.param.Pt_snk ) 
      {
	string source_sink = Mass_12 + ".D" + src_wvf_param_s[i2] + ".P"  + ".SP";
	int nn;

	files.had[i].current_SP.vector.resize(spec.had[i].current_SP.vector.size());
	nn = spec.had[i].current_SP.vector.size();
	for(int n = 0; n < nn; ++n)
	  files.had[i].current_SP.vector[n].filename = vector_particle[n] + source_sink;

	files.had[i].current_SP.axial.resize(spec.had[i].current_SP.axial.size());
	nn = spec.had[i].current_SP.axial.size();
	for(int n = 0; n < nn; ++n)
	  files.had[i].current_SP.axial[n].filename = axial_particle[n] + source_sink;
      }

      if ( spec.had[i].Wl_src && spec.param.Pt_snk ) 
      {
	string source_sink = Mass_12 + ".W.P.WP";
	int nn;

	files.had[i].current_WP.vector.resize(spec.had[i].current_WP.vector.size());
	nn = spec.had[i].current_WP.vector.size();
	for(int n = 0; n < nn; ++n)
	  files.had[i].current_WP.vector[n].filename = vector_particle[n] + source_sink;

	files.had[i].current_WP.axial.resize(spec.had[i].current_WP.axial.size());
	nn = spec.had[i].current_WP.axial.size();
	for(int n = 0; n < nn; ++n)
	  files.had[i].current_WP.axial[n].filename = axial_particle[n] + source_sink;
      }
    }  // end for i
  }



  //
  // Construct hybrid meson names
  //
  if (spec.param.HybMesP)
  {
    Array<string> hybmes_particle(15);

    hybmes_particle[0]  = "hpion";
    hybmes_particle[1]  = "hrho_x";
    hybmes_particle[2]  = "hrho_y";
    hybmes_particle[3]  = "hrho_z";
    hybmes_particle[4]  = "ex0pm";
    hybmes_particle[5]  = "ex0mm";
    hybmes_particle[6]  = "ex1mp_x";
    hybmes_particle[7]  = "ex1mp_y";
    hybmes_particle[8]  = "ex1mp_z";
    hybmes_particle[9]  = "ex1mp_x";
    hybmes_particle[10] = "ex1mp_y";
    hybmes_particle[11] = "ex1mp_z";
    hybmes_particle[12] = "ex1mp_x";
    hybmes_particle[13] = "ex1mp_y";
    hybmes_particle[14] = "ex1mp_z";

    Array<string> hybmes_smear_state(15);

    hybmes_smear_state[0]  = "_1" ;
    hybmes_smear_state[1]  = "_1" ;
    hybmes_smear_state[2]  = "_1" ;
    hybmes_smear_state[3]  = "_1" ;
    hybmes_smear_state[4]  = "_1" ;
    hybmes_smear_state[5]  = "_1" ;
    hybmes_smear_state[6]  = "_1" ;
    hybmes_smear_state[7]  = "_1" ;
    hybmes_smear_state[8]  = "_1" ;
    hybmes_smear_state[9]  = "_2" ;
    hybmes_smear_state[10] = "_2" ;
    hybmes_smear_state[11] = "_2" ;
    hybmes_smear_state[12] = "_3" ;
    hybmes_smear_state[13] = "_3" ;
    hybmes_smear_state[14] = "_3" ;

    /*
     * Create hybrid meson source-sink smearing names
     */
    cout << "Creating hybrid meson file names" << endl;


    // For now assume only diagonal valence mesons
    for(int i = 0; i < spec.had.size(); ++i) 
    {
      int i1=i, i2=i;
      bool diag = true;
      string Mass_12 = ".D" + Mass_s[i1];
      string wvf_param12_s = "";   // ignore for the moment, should fill in

      if ( spec.had[i].Pt_src && spec.param.Pt_snk ) 
      {
	files.had[i].hybmes_PP.resize(spec.had[i].hybmes_PP.size());

	for(int m=0; m < spec.had[i].hybmes_PP.size(); ++m)
	{
	  string ms = hybmes_smear_state[m];
	
	  files.had[i].hybmes_PP[m].mom.resize(spec.had[i].hybmes_PP[m].mom.size());
	
	  for(int n = 0; n < spec.had[i].hybmes_PP[m].mom.size(); ++n)
	  {
	    files.had[i].hybmes_PP[m].mom[n].filename = hybmes_particle[m] 
	      + mom_name(spec.had[i].hybmes_PP[m].mom[n].sink_mom) 
	      + Mass_12 + ".P" + ms + ".P" + ms + ".PP";
	  }
	}
      }

      if ( spec.had[i].Pt_src && spec.param.Sl_snk ) 
      {
	files.had[i].hybmes_PS.resize(spec.had[i].hybmes_PS.size());

	for(int m=0; m < spec.had[i].hybmes_PS.size(); ++m)
	{
	  string ms = hybmes_smear_state[m];

	  files.had[i].hybmes_PS[m].mom.resize(spec.had[i].hybmes_PS[m].mom.size());

	  for(int n = 0; n < spec.had[i].hybmes_PS[m].mom.size(); ++n)
	  {
	    files.had[i].hybmes_PS[m].mom[n].filename = hybmes_particle[m] 
	      + mom_name(spec.had[i].hybmes_PS[m].mom[n].sink_mom) 
	      + Mass_12 + ".P" + ms + ".D" + snk_wvf_param_s[i2] + ms + ".PS";
	  }
	}
      }

      if ( spec.had[i].Pt_src && spec.param.Wl_snk ) 
      {
	files.had[i].hybmes_PW.resize(spec.had[i].hybmes_PW.size());

	for(int m=0; m < spec.had[i].hybmes_PW.size(); ++m)
	{
	  string ms = hybmes_smear_state[m];
	
	  files.had[i].hybmes_PW[m].mom.resize(spec.had[i].hybmes_PW[m].mom.size());
	
	  for(int n = 0; n < spec.had[i].hybmes_PW[m].mom.size(); ++n)
	  {
	    files.had[i].hybmes_PW[m].mom[n].filename = hybmes_particle[m] 
	      + mom_name(spec.had[i].hybmes_PW[m].mom[n].sink_mom) 
	      + Mass_12 + ".P" + ms + ".W" + ms + ".PW";
	  }
	}
      }

      if ( spec.had[i].Sl_src && spec.param.Pt_snk ) 
      {
	files.had[i].hybmes_SP.resize(spec.had[i].hybmes_SP.size());

	for(int m=0; m < spec.had[i].hybmes_SP.size(); ++m)
	{
	  string ms = hybmes_smear_state[m];

	  files.had[i].hybmes_SP[m].mom.resize(spec.had[i].hybmes_SP[m].mom.size());

	  for(int n = 0; n < spec.had[i].hybmes_SP[m].mom.size(); ++n)
	  {
	    files.had[i].hybmes_SP[m].mom[n].filename = hybmes_particle[m] 
	      + mom_name(spec.had[i].hybmes_SP[m].mom[n].sink_mom) 
	      + Mass_12 + ".D" + src_wvf_param_s[i2] + ms + ".P" + ms + ".SP";

//	    cout << spec.had[i].hybmes_SP[m].mom[n].filename << endl;
	  }
	}
      }

      if ( spec.had[i].Sl_src && spec.param.Sl_snk ) 
      {
	files.had[i].hybmes_SS.resize(spec.had[i].hybmes_SS.size());

	for(int m=0; m < spec.had[i].hybmes_SS.size(); ++m)
	{
	  string ms = hybmes_smear_state[m];

	  files.had[i].hybmes_SS[m].mom.resize(spec.had[i].hybmes_SS[m].mom.size());

	  for(int n = 0; n < spec.had[i].hybmes_SS[m].mom.size(); ++n)
	  {
	    files.had[i].hybmes_SS[m].mom[n].filename = hybmes_particle[m] 
	      + mom_name(spec.had[i].hybmes_SS[m].mom[n].sink_mom) 
	      + Mass_12 + ".D" + src_wvf_param_s[i1] + ms + ".D" + snk_wvf_param_s[i2] + ms + ".SS";

//	    cout << spec.had[i].hybmes_SP[m].mom[n].filename << endl;
	  }
	}
      }

      if ( spec.had[i].Sl_src && spec.param.Wl_snk ) 
      {
	files.had[i].hybmes_SW.resize(spec.had[i].hybmes_SW.size());

	for(int m=0; m < spec.had[i].hybmes_SW.size(); ++m)
	{
	  string ms = hybmes_smear_state[m];

	  files.had[i].hybmes_SW[m].mom.resize(spec.had[i].hybmes_SW[m].mom.size());

	  for(int n = 0; n < spec.had[i].hybmes_SW[m].mom.size(); ++n)
	  {
	    files.had[i].hybmes_SW[m].mom[n].filename = hybmes_particle[m] 
	      + mom_name(spec.had[i].hybmes_SW[m].mom[n].sink_mom) 
	      + Mass_12 + ".D" + src_wvf_param_s[i2] + ms + ".W" + ms + ".SW";

//	    cout << spec.had[i].hybmes_SW[m].mom[n].filename << endl;
	  }
	}
      }

      if ( spec.had[i].Wl_src && spec.param.Pt_snk ) 
      {
	files.had[i].hybmes_WP.resize(spec.had[i].hybmes_WP.size());

	for(int m=0; m < spec.had[i].hybmes_WP.size(); ++m)
	{
	  string ms = hybmes_smear_state[m];
	
	  files.had[i].hybmes_WP[m].mom.resize(spec.had[i].hybmes_WP[m].mom.size());
	
	  for(int n = 0; n < spec.had[i].hybmes_WP[m].mom.size(); ++n)
	  {
	    files.had[i].hybmes_WP[m].mom[n].filename = hybmes_particle[m] 
	      + mom_name(spec.had[i].hybmes_WP[m].mom[n].sink_mom) 
	      + Mass_12 + ".W" + ms + ".P" + ms + ".WP";
	  }
	}
      }

      if ( spec.had[i].Wl_src && spec.param.Sl_snk ) 
      {
	files.had[i].hybmes_WS.resize(spec.had[i].hybmes_WS.size());

	for(int m=0; m < spec.had[i].hybmes_WS.size(); ++m)
	{
	  string ms = hybmes_smear_state[m];

	  files.had[i].hybmes_WS[m].mom.resize(spec.had[i].hybmes_WS[m].mom.size());

	  for(int n = 0; n < spec.had[i].hybmes_WS[m].mom.size(); ++n)
	  {
	    files.had[i].hybmes_WS[m].mom[n].filename = hybmes_particle[m] 
	      + mom_name(spec.had[i].hybmes_WS[m].mom[n].sink_mom) 
	      + Mass_12 + ".W" + ms + ".D" + snk_wvf_param_s[i2] + ms + ".WS";
	  }
	}
      }

      if ( spec.had[i].Wl_src && spec.param.Wl_snk ) 
      {
	files.had[i].hybmes_WW.resize(spec.had[i].hybmes_WW.size());

	for(int m=0; m < spec.had[i].hybmes_WW.size(); ++m)
	{
	  string ms = hybmes_smear_state[m];
	
	  files.had[i].hybmes_WW[m].mom.resize(spec.had[i].hybmes_WW[m].mom.size());
	
	  for(int n = 0; n < spec.had[i].hybmes_WW[m].mom.size(); ++n)
	  {
	    files.had[i].hybmes_WW[m].mom[n].filename = hybmes_particle[m] 
	      + mom_name(spec.had[i].hybmes_WW[m].mom[n].sink_mom) 
	      + Mass_12 + ".W" + ms + ".W" + ms + ".WW";
	  }
	}
      }
    } 
  }


  //
  // Construct baryon names
  //
  if (spec.param.BaryonP)
  {
    Array<string> baryon_particle(22);

    baryon_particle[0]  = "proton";
    baryon_particle[1]  = "lambda";
    baryon_particle[2]  = "delta";
    baryon_particle[3]  = "proton";
    baryon_particle[4]  = "lambda";
    baryon_particle[5]  = "delta";
    baryon_particle[6]  = "proton";
    baryon_particle[7]  = "lambda";
    baryon_particle[8]  = "delta";
    baryon_particle[9]  = "proton";
    baryon_particle[10] = "proton";
    baryon_particle[11] = "proton";
    baryon_particle[12]  = "delta_x";
    baryon_particle[13]  = "delta_y";
    baryon_particle[14]  = "delta_z";
    baryon_particle[15]  = "delta_x";
    baryon_particle[16]  = "delta_y";
    baryon_particle[17]  = "delta_z";
    baryon_particle[18]  = "delta_x";
    baryon_particle[19]  = "delta_y";
    baryon_particle[20]  = "delta_z";
    baryon_particle[21]  = "proton_negpar";

    Array<string> baryon_smear_state(22);

    baryon_smear_state[0]  = "_1" ;
    baryon_smear_state[1]  = "_1" ;
    baryon_smear_state[2]  = "_1" ;
    baryon_smear_state[3]  = "_2" ;
    baryon_smear_state[4]  = "_2" ;
    baryon_smear_state[5]  = "_2" ;
    baryon_smear_state[6]  = "_3" ;
    baryon_smear_state[7]  = "_3" ;
    baryon_smear_state[8]  = "_3" ;
    baryon_smear_state[9]  = "_4" ;
    baryon_smear_state[10] = "_5" ;
    baryon_smear_state[11] = "_6" ;
    baryon_smear_state[12]  = "_4" ;
    baryon_smear_state[13]  = "_4" ;
    baryon_smear_state[14]  = "_4" ;
    baryon_smear_state[15] = "_5" ;
    baryon_smear_state[16] = "_5" ;
    baryon_smear_state[17] = "_5" ;
    baryon_smear_state[18] = "_6" ;
    baryon_smear_state[19] = "_6" ;
    baryon_smear_state[20] = "_6" ;
    baryon_smear_state[21] = "_3" ;

    /*
     * Create Baryon source-sink smearing names
     */
    cout << "Creating baryon file names" << endl;

    // For now assume only diagonal valence baryons
    for(int i = 0; i < spec.had.size(); ++i) 
    {
      int i1=i, i2=i, i3=i;
      bool diag = true;
      string Mass_123 = ".D" + Mass_s[i1];
      string wvf_param123_s = "";   // ignore for the moment, should fill in

      if ( spec.had[i].Pt_src && spec.param.Pt_snk ) 
      {
//	files.had[i].baryon_PP.resize(spec.had[i].baryon_PP.size());
	files.had[i].baryon_PP.resize(baryon_particle.size());

	for(int m=0; m < files.had[i].baryon_PP.size(); ++m)
	{
	  string ms = baryon_smear_state[m];
	
	  files.had[i].baryon_PP[m].mom.resize(spec.had[i].baryon_PP[0].mom.size());
	
	  for(int n = 0; n < spec.had[i].baryon_PP[0].mom.size(); ++n)
	  {
	    files.had[i].baryon_PP[m].mom[n].filename = baryon_particle[m] 
	      + mom_name(spec.had[i].baryon_PP[0].mom[n].sink_mom) 
	      + Mass_123 + ".P" + ms + ".P" + ms + ".PP";
	  }
	}
      }

      if ( diag )
      {
	if ( spec.had[i].Pt_src && spec.param.Sl_snk ) 
	{
//	  files.had[i].baryon_PS.resize(spec.had[i].baryon_PS.size());
	  files.had[i].baryon_PS.resize(baryon_particle.size());

	  for(int m=0; m < files.had[i].baryon_PS.size(); ++m)
	  {
	    string ms = baryon_smear_state[m];

	    files.had[i].baryon_PS[m].mom.resize(spec.had[i].baryon_PS[0].mom.size());

	    for(int n = 0; n < spec.had[i].baryon_PS[0].mom.size(); ++n)
	    {
	      files.had[i].baryon_PS[m].mom[n].filename = baryon_particle[m] 
		+ mom_name(spec.had[i].baryon_PS[0].mom[n].sink_mom) 
		+ Mass_123 + ".P" + ms + ".D" + snk_wvf_param_s[i2] + ms + ".PS";
	    }
	  }
	}

	if ( spec.had[i].Pt_src && spec.param.Wl_snk ) 
	{
//	  files.had[i].baryon_PW.resize(spec.had[i].baryon_PW.size());
	  files.had[i].baryon_PW.resize(baryon_particle.size());

	  for(int m=0; m < files.had[i].baryon_PW.size(); ++m)
	  {
	    string ms = baryon_smear_state[m];
	
	    files.had[i].baryon_PW[m].mom.resize(spec.had[i].baryon_PW[0].mom.size());
	
	    for(int n = 0; n < spec.had[i].baryon_PW[0].mom.size(); ++n)
	    {
	      files.had[i].baryon_PW[m].mom[n].filename = baryon_particle[m] 
		+ mom_name(spec.had[i].baryon_PW[0].mom[n].sink_mom) 
		+ Mass_123 + ".P" + ms + ".W" + ms + ".PW";
	    }
	  }
	}
	
	if ( spec.had[i].Sl_src && spec.param.Pt_snk ) 
	{
//	  files.had[i].baryon_SP.resize(spec.had[i].baryon_SP.size());
	  files.had[i].baryon_SP.resize(baryon_particle.size());

	  for(int m=0; m < files.had[i].baryon_SP.size(); ++m)
	  {
	    string ms = baryon_smear_state[m];

	    files.had[i].baryon_SP[m].mom.resize(spec.had[i].baryon_SP[0].mom.size());

	    for(int n = 0; n < spec.had[i].baryon_SP[0].mom.size(); ++n)
	    {
	      files.had[i].baryon_SP[m].mom[n].filename = baryon_particle[m] 
		+ mom_name(spec.had[i].baryon_SP[0].mom[n].sink_mom) 
		+ Mass_123 + ".D" + src_wvf_param_s[i2] + ms + ".P" + ms + ".SP";
	    }
	  }
	}

	if ( spec.had[i].Sl_src && spec.param.Sl_snk ) 
	{
//	  files.had[i].baryon_SS.resize(spec.had[i].baryon_SS.size());
	  files.had[i].baryon_SS.resize(baryon_particle.size());

	  for(int m=0; m < files.had[i].baryon_SS.size(); ++m)
	  {
	    string ms = baryon_smear_state[m];

	    files.had[i].baryon_SS[m].mom.resize(spec.had[i].baryon_SS[0].mom.size());

	    for(int n = 0; n < spec.had[i].baryon_SS[0].mom.size(); ++n)
	    {
	      files.had[i].baryon_SS[m].mom[n].filename = baryon_particle[m] 
		+ mom_name(spec.had[i].baryon_SS[0].mom[n].sink_mom) 
		+ Mass_123 + ".D" + src_wvf_param_s[i1] + ms + ".D" + snk_wvf_param_s[i2] + ms + ".SS";
	    }
	  }
	}

	if ( spec.had[i].Sl_src && spec.param.Wl_snk ) 
	{
//	  files.had[i].baryon_SW.resize(spec.had[i].baryon_SW.size());
	  files.had[i].baryon_SW.resize(baryon_particle.size());

	  for(int m=0; m < files.had[i].baryon_SW.size(); ++m)
	  {
	    string ms = baryon_smear_state[m];

	    files.had[i].baryon_SW[m].mom.resize(spec.had[i].baryon_SW[0].mom.size());

	    for(int n = 0; n < spec.had[i].baryon_SW[0].mom.size(); ++n)
	    {
	      files.had[i].baryon_SW[m].mom[n].filename = baryon_particle[m] 
		+ mom_name(spec.had[i].baryon_SW[0].mom[n].sink_mom) 
		+ Mass_123 + ".D" + src_wvf_param_s[i2] + ms + ".W" + ms + ".SW";
	    }
	  }
	}

	if ( spec.had[i].Wl_src && spec.param.Pt_snk ) 
	{
//	  files.had[i].baryon_WP.resize(spec.had[i].baryon_WP.size());
	  files.had[i].baryon_WP.resize(baryon_particle.size());

	  for(int m=0; m < files.had[i].baryon_WP.size(); ++m)
	  {
	    string ms = baryon_smear_state[m];
	
	    files.had[i].baryon_WP[m].mom.resize(spec.had[i].baryon_WP[0].mom.size());
	
	    for(int n = 0; n < spec.had[i].baryon_WP[0].mom.size(); ++n)
	    {
	      files.had[i].baryon_WP[m].mom[n].filename = baryon_particle[m] 
		+ mom_name(spec.had[i].baryon_WP[0].mom[n].sink_mom) 
		+ Mass_123 + ".W" + ms + ".P" + ms + ".WP";
	    }
	  }
	}

	if ( spec.had[i].Wl_src && spec.param.Sl_snk ) 
	{
//	  files.had[i].baryon_WS.resize(spec.had[i].baryon_WS.size());
	  files.had[i].baryon_WS.resize(baryon_particle.size());

	  for(int m=0; m < files.had[i].baryon_WS.size(); ++m)
	  {
	    string ms = baryon_smear_state[m];

	    files.had[i].baryon_WS[m].mom.resize(spec.had[i].baryon_WS[0].mom.size());

	    for(int n = 0; n < spec.had[i].baryon_WS[0].mom.size(); ++n)
	    {
	      files.had[i].baryon_WS[m].mom[n].filename = baryon_particle[m] 
		+ mom_name(spec.had[i].baryon_WS[0].mom[n].sink_mom) 
		+ Mass_123 + ".W" + ms + ".D" + snk_wvf_param_s[i2] + ms + ".WS";
	    }
	  }
	}

	if ( spec.had[i].Wl_src && spec.param.Wl_snk ) 
	{
//	  files.had[i].baryon_WW.resize(spec.had[i].baryon_WW.size());
	  files.had[i].baryon_WW.resize(baryon_particle.size());

	  for(int m=0; m < files.had[i].baryon_WW.size(); ++m)
	  {
	    string ms = baryon_smear_state[m];
	
	    files.had[i].baryon_WW[m].mom.resize(spec.had[i].baryon_WW[0].mom.size());
	
	    for(int n = 0; n < spec.had[i].baryon_WW[0].mom.size(); ++n)
	    {
	      files.had[i].baryon_WW[m].mom[n].filename = baryon_particle[m] 
		+ mom_name(spec.had[i].baryon_WW[0].mom[n].sink_mom) 
		+ Mass_123 + ".W" + ms + ".W" + ms + ".WW";
	    }
	  }
	}

      } 
      else    // if (! diag )   NOTE: the code will never get to this part for now
      {
	if ( spec.had[i].Pt_src && spec.param.Sl_snk ) 
	{
	  files.had[i].baryon_PS.resize(spec.had[i].baryon_PS.size());

	  for(int m=0; m < files.had[i].baryon_PS.size(); ++m)
	  {
	    string ms = baryon_smear_state[m];

	    files.had[i].baryon_PS[m].mom.resize(spec.had[i].baryon_PS[0].mom.size());

	    for(int n = 0; n < spec.had[i].baryon_PS[0].mom.size(); ++n)
	    {
	      files.had[i].baryon_PS[m].mom[n].filename = baryon_particle[m] 
		+ mom_name(spec.had[i].baryon_PS[0].mom[n].sink_mom) 
		+ Mass_123 + ".P" + ms + wvf_param123_s + ms + ".PS";
	    }
	  }
	}

	if ( spec.had[i].Sl_src && spec.param.Pt_snk ) 
	{
	  files.had[i].baryon_SP.resize(spec.had[i].baryon_SP.size());

	  for(int m=0; m < files.had[i].baryon_SP.size(); ++m)
	  {
	    string ms = baryon_smear_state[m];

	    files.had[i].baryon_SP[m].mom.resize(spec.had[i].baryon_SP[0].mom.size());

	    for(int n = 0; n < spec.had[i].baryon_SP[0].mom.size(); ++n)
	    {
	      files.had[i].baryon_SP[m].mom[n].filename = baryon_particle[m] 
		+ mom_name(spec.had[i].baryon_SP[0].mom[n].sink_mom) 
		+ Mass_123 + wvf_param123_s + ms + ".P" + ms + ".SP";
	    }
	  }
	}

	if ( spec.had[i].Sl_src && spec.param.Sl_snk ) 
	{
	  files.had[i].baryon_SS.resize(spec.had[i].baryon_SS.size());

	  for(int m=0; m < files.had[i].baryon_SS.size(); ++m)
	  {
	    string ms = baryon_smear_state[m];

	    files.had[i].baryon_SS[m].mom.resize(spec.had[i].baryon_SS[0].mom.size());

	    for(int n = 0; n < spec.had[i].baryon_SS[0].mom.size(); ++n)
	    {
	      files.had[i].baryon_SS[m].mom[n].filename = baryon_particle[m] 
		+ mom_name(spec.had[i].baryon_SS[0].mom[n].sink_mom) 
		+ Mass_123 + wvf_param123_s + ms + wvf_param123_s + ms + ".SS";
	    }
	  }
	}

      }  // end if (diag == 1)
    }  // end for i
  }
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

  string xml_in_root = "/spectrum_w";


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
  cout << "spec.had.meson_PP.size = " << spec.had[0].meson_PP.size() << endl;
  cout << "spec.had.meson_PS.size = " << spec.had[0].meson_PS.size() << endl;
  cout << "spec.had.meson_SP.size = " << spec.had[0].meson_SP.size() << endl;
  cout << "spec.had.meson_SS.size = " << spec.had[0].meson_SS.size() << endl;
//  cout << "spec.had.meson_SP.mom.size = " << spec.had[0].meson_SP[0].mom.size() << endl;
//  cout << "spec.had.meson_SP.mom.mesprop.size = " << spec.had[0].meson_SP[0].mom[0].mesprop.size() << endl;
  cout << "spec.had.hybmes_SP.size = " << spec.had[0].hybmes_SP.size() << endl;
  
  if (spec.param.wvf_param.size() != spec.had.size())
  {
    cerr << "wvf_param and Mass array do not agree in size";
    exit(1);
  }


  // We care about the output version
  switch (spec.output_version.out_version)
  {
  case 9:
  case 10:
  case 11:
  case 12:
  case 13:
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

  exit(0);
}
