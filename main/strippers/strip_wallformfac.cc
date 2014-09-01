// $Id: strip_wallformfac.cc,v 2.0 2008/12/05 04:44:04 edwards Exp $

#include <sys/time.h>   // for timings

#include "strippers.h"

using namespace Strippers;
using namespace std;

const int Nd = 4;

/*
 * Input 
 */
struct Output_version_t
{
  int out_version;
};

struct Param_t
{
  int          mom2_max;
  bool         wall_source;         // use wall source or wall sink
  Array<int>   nrow;
};

/*
 * Structures for hadron parts
 */
struct WallFormFac_momenta_t
{
  int              inser_mom_num;
  Array<int>       inser_mom;
  Array<Complex>   local_current;
  Array<Complex>   nonlocal_current;
};

struct WallFormFac_insertion_t
{
  int              gamma_ctr;
  int              mu;
  int              gamma_value;
  Array<WallFormFac_momenta_t> momenta;
};

struct WallFormFac_projector_t
{
  int              proj_ctr;
  string           proj_name;
  Array<WallFormFac_insertion_t>  insertion;
};

struct WallFormFac_lorentz_t
{
  int              lorentz_ctr;
  int              snk_gamma;
  int              src_gamma;
  Array<WallFormFac_projector_t>  projector;
};

struct WallFormFac_formfac_t
{
  int              formfac_ctr;
  string           formfac_name;
  Array<WallFormFac_lorentz_t>  lorentz;
};

struct WallFormFac_quark_t
{
  int              quark_ctr;
  string           quark_name;
  Array<WallFormFac_formfac_t>  formfac;
};

struct WallFormFac_formfacs_t
{
  string           subroutine;
  Array<WallFormFac_quark_t>  quark;
};

struct WallFormFac_bar_t
{
  int                     formfac_value;
  string                  formfac_type;
  WallFormFac_formfacs_t  formfacs;
};


struct WallFormFac_t
{
  int          out_version;
  int          mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.
  bool         wall_source;         // use wall source or wall sink
  Array<int> nrow;
  Array<WallFormFac_bar_t> bar;
};


/*
 * Similar structure to Wallformfac_t that holds filenames
 * NOTE: this cannot be merged into Wallformfac_t since that
 * is volatile - upon reading this struct old data is
 * wiped out
 */
struct File_t
{
  string    filename;   // used by file writer below
};

struct File_momenta_t
{
  File_t local_current;
  File_t nonlocal_current;
};

struct File_insertion_t
{
  Array<File_momenta_t> momenta;
};

struct File_projector_t
{
  Array<File_insertion_t>  insertion;
};

struct File_lorentz_t
{
  Array<File_projector_t>  projector;
};

struct File_formfac_t
{
  Array<File_lorentz_t>  lorentz;
};

struct File_quark_t
{
  Array<File_formfac_t>  formfac;
};

struct File_formfacs_t
{
  Array<File_quark_t>  quark;
};

struct File_bar_t
{
  File_formfacs_t  formfacs;
};

struct File_wallformfac_t
{
  int nprop;
  Array<File_bar_t>  bar;
};


// Read a complex
void read(XMLReader& xml, const string& path, Complex& com)
{
  XMLReader complextop(xml, path);
  read(complextop, "re", com.re);
  read(complextop, "im", com.im);
}


// Read formfactors
void read(BinaryReader& bin, Complex& p)
{
  read(bin, p.re);
  read(bin, p.im);
}



/*
 * Read the parameters
 */
void read(XMLReader& xml, const string& path, Output_version_t& input)
{
  XMLReader inputtop(xml, path);

  // Not setting context here as it is only the one read.
  read(inputtop, "out_version", input.out_version);
}


//! Parameter input
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader paramtop(xml, path);

  int version;
  read(paramtop, "version", version);

  switch (version) 
  {
    /**************************************************************************/
  case 2:
    param.wall_source = false;
    break;

    /**************************************************************************/
  case 3:
    read(paramtop, "wall_source", param.wall_source);
    break;

  default:
    /**************************************************************************/
    cerr << "Input parameter version " << version 
		<< " unsupported." << endl;
    throw;
  }

  read(paramtop, "mom2_max", param.mom2_max);
  read(paramtop, "nrow", param.nrow);
}



// Read a momenta struct
void read(XMLReader& xml, const string& path, WallFormFac_momenta_t& mom)
{
  XMLReader top(xml, path);

  read(top, "inser_mom", mom.inser_mom);

  read(top, "local_cur3ptfn", mom.local_current);

  // Check for existence first
  string xpath ="nonlocal_cur3ptfn";
  if (top.count(xpath) != 0) {
    read(top, xpath, mom.nonlocal_current);
  }
}

//! Wallformfac insertion writer
void read(XMLReader& xml, const string& path, WallFormFac_insertion_t& header)
{
  XMLReader top(xml, path);

  read(top, "gamma_ctr", header.gamma_ctr);
  read(top, "mu", header.mu);
  read(top, "gamma_value", header.gamma_value);
  read(top, "Momenta", header.momenta);
}

//! Wallformfac projector writer
void read(XMLReader& xml, const string& path, WallFormFac_projector_t& header)
{
  XMLReader top(xml, path);

  read(top, "proj_ctr", header.proj_ctr);
  read(top, "proj_name", header.proj_name);
  read(top, "Insertion", header.insertion);
}

//! Wallformfac lorentz writer
void read(XMLReader& xml, const string& path, WallFormFac_lorentz_t& header)
{
  XMLReader top(xml, path);

  read(top, "lorentz_ctr", header.lorentz_ctr);
  read(top, "snk_gamma", header.snk_gamma);
  read(top, "src_gamma", header.src_gamma);
  read(top, "Projector", header.projector);
}

//! Wallformfac formfac writer
void read(XMLReader& xml, const string& path, WallFormFac_formfac_t& header)
{
  XMLReader top(xml, path);

  read(top, "formfac_ctr", header.formfac_ctr);
  read(top, "formfac_name", header.formfac_name);
  read(top, "Lorentz", header.lorentz);
}

//! Wallformfac quark writer
void read(XMLReader& xml, const string& path, WallFormFac_quark_t& header)
{
  XMLReader top(xml, path);

  read(top, "quark_ctr", header.quark_ctr);
  read(top, "quark_name", header.quark_name);
  read(top, "FormFac", header.formfac);
}

//! WallFormFac writer
void read(XMLReader& xml, const string& path, WallFormFac_formfacs_t& header)
{
  XMLReader top(xml, path);

  read(top, "subroutine", header.subroutine);
  read(top, "Quark", header.quark);
}

void read(XMLReader& xml, const string& path, WallFormFac_bar_t& header)
{
  XMLReader top(xml, path);

  read(top, "formfac_type", header.formfac_type);
  read(top, "formfac_value", header.formfac_value);
  read(top, "WallFormFac", header.formfacs);
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, WallFormFac_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // Parameters for source construction
  // Not setting context here as it is only the one read.
    Output_version_t output_version;
    read(inputtop, "Output_version", output_version);
    input.out_version = output_version.out_version;

    // Parameters for source construction
    Param_t param;
    read(inputtop, "Input/WallFormFac/Param", param);
    input.mom2_max = param.mom2_max;
    input.wall_source = param.wall_source;
    input.nrow = param.nrow;

    // Read the data
    read(inputtop, "Wilson_3Pt_fn_measurements", input.bar);
  }
  catch(const string& e)
  {
    cerr << "Error reading in wallformfac: " << e << endl;
    throw;
  }
}



//! Wallformfac momenta writer
void read(BinaryReader& bin, WallFormFac_momenta_t& header)
{
  int magic;
  read(bin, magic);
  if (magic != 20301)
  {
    cerr << "read(Momenta_t): magic number invalid" << endl;
    throw;
  }
  read(bin, header.inser_mom_num);
  read(bin, header.inser_mom);
  read(bin, header.local_current);
//  cout << "local_current.size() = " << header.local_current.size() << endl;
  read(bin, header.nonlocal_current);
//  cout << "nonlocal_current.size() = " << header.nonlocal_current.size() << endl;
}

//! Wallformfac insertion writer
void read(BinaryReader& bin, WallFormFac_insertion_t& header)
{
  read(bin, header.gamma_ctr);
  read(bin, header.mu);
  read(bin, header.gamma_value);
  read(bin, header.momenta);
}

//! Wallformfac projector writer
void read(BinaryReader& bin, WallFormFac_projector_t& header)
{
  read(bin, header.proj_ctr);
  read(bin, header.proj_name, 100);
//  cout << "proj_name= " << header.proj_name << endl;
  read(bin, header.insertion);
}

//! Wallformfac lorentz writer
void read(BinaryReader& bin, WallFormFac_lorentz_t& header)
{
  read(bin, header.lorentz_ctr);
  read(bin, header.snk_gamma);
  read(bin, header.src_gamma);
  read(bin, header.projector);
}

//! Wallformfac formfac writer
void read(BinaryReader& bin, WallFormFac_formfac_t& header)
{
  read(bin, header.formfac_ctr);
  read(bin, header.formfac_name, 100);
//  cout << "formfac_name= " << header.formfac_name << endl;
  read(bin, header.lorentz);
}

//! Wallformfac quark writer
void read(BinaryReader& bin, WallFormFac_quark_t& header)
{
  read(bin, header.quark_ctr);
  read(bin, header.quark_name, 100);
//  cout << "quark_name= " << header.quark_name << endl;
  read(bin, header.formfac);
}

//! WallFormFac writer
void read(BinaryReader& bin, WallFormFac_formfacs_t& header)
{
  read(bin, header.subroutine, 100);
//  cout << "subroutine_name= " << header.subroutine << endl;
  read(bin, header.quark);
}

//! WallFormFac writer
void read(BinaryReader& bin, WallFormFac_bar_t& header)
{
  read(bin, header.formfac_value);
  read(bin, header.formfac_type, 100);
//  cout << "formfac_type= " << header.formfac_type << endl;
  read(bin, header.formfacs);
}

//! WallFormFac writer
void read(BinaryReader& bin, WallFormFac_t& header)
{
  try
  {
    read(bin, header.out_version);
    read(bin, header.mom2_max);

    switch (header.out_version) 
    {
      /**************************************************************************/
    case 3:
      header.wall_source = false;
      break;

      /**************************************************************************/
    case 4:
      read(bin, header.wall_source);
      break;

    default:
      /**************************************************************************/
      cerr << "Output parameter version " << header.out_version 
	   << " unsupported." << endl;
      throw;
    }

    read(bin, header.nrow);
    read(bin, header.bar);
  }
  catch(const string& e)
  {
    cerr << "Error reading binary in wallformfac: " << e << endl;
    throw;
  }
}



// Write formfactors
ostream& operator<<(ostream& s, const Array<Complex>& p)
{
  for(int i=0; i < p.size(); ++i)
    s << i << " " << p[i].re << " " << p[i].im << "\n";
  return s;
}

// Write stuff
ostream& operator<<(ostream& s, const Array<int>& p)
{
  for(int i=0; i < p.size(); ++i)
    s << " " << p[i];
  s << endl;
  return s;
}


string mom_name(const Array<int>& sink_mom)
{
  ostringstream s;

  // Always print momenta structure
  s << "_qx" << sink_mom[0] << "_qy" << sink_mom[1] << "_qz" << sink_mom[2];

  return s.str();
}


void print_file(const string& filename, 
		const Array<Complex>& current, 
		int nprop, 
		const Array<int>& nrow, 
		bool first)
{
  if (first)
  {
    ofstream file(filename.c_str());

    // Write header line
// FIX THIS - RATS - need to GET j_decay BACK IN HERE!!
    file << nprop << " " << nrow[Nd-1] << " 1 " 
	 << nrow[0] << " 2" << endl;

    file << current;   // Write the file
    file.close();
  }
  else
  {
    ofstream file(filename.c_str(), std::ios_base::app);
    file << current;
    file.close();
  }
}


void print_files(const string& formfac_type,
		 const File_wallformfac_t& files, const WallFormFac_t& bar3pt, bool first)
{
  for(int b=0; b < files.bar.size(); ++b)
  {
    const WallFormFac_bar_t& bar = bar3pt.bar[b];
    const File_bar_t& file_bar = files.bar[b];

    if ((formfac_type != bar.formfac_type) && (formfac_type != "all"))
      continue;

    const WallFormFac_formfacs_t& formfacs = bar.formfacs;
    const File_formfacs_t& file_formfacs = file_bar.formfacs;

    for(int q=0; q < file_formfacs.quark.size(); ++q)
    {
      const WallFormFac_quark_t& quark = formfacs.quark[q];
      const File_quark_t& file_quark = file_formfacs.quark[q];

      for(int f=0; f < file_quark.formfac.size(); ++f)
      {
	const WallFormFac_formfac_t& formfac = quark.formfac[f];
	const File_formfac_t& file_formfac = file_quark.formfac[f];
	
	for(int l=0; l < file_formfac.lorentz.size(); ++l)
	{
	  const WallFormFac_lorentz_t& lorentz = formfac.lorentz[l];
	  const File_lorentz_t& file_lorentz = file_formfac.lorentz[l];

	  for(int p=0; p < file_lorentz.projector.size(); ++p)
	  {
	    const WallFormFac_projector_t& projector = lorentz.projector[p];
	    const File_projector_t& file_projector = file_lorentz.projector[p];

	    for(int g=0; g < file_projector.insertion.size(); ++g)
	    {
	      const WallFormFac_insertion_t& insertion = projector.insertion[g];
	      const File_insertion_t& file_insertion = file_projector.insertion[g];

	      for(int n = 0; n < file_insertion.momenta.size(); ++n)
	      {
		const WallFormFac_momenta_t& momenta = insertion.momenta[n];
		const File_momenta_t& file_momenta = file_insertion.momenta[n];

		// The local current should always be there, but check anyway
		if (momenta.local_current.size() > 0)
		  print_file(file_momenta.local_current.filename,
			     momenta.local_current,
			     files.nprop, bar3pt.nrow, first);
		
		// The non-local current may not always be present, so check
		if (momenta.nonlocal_current.size() > 0)
		  print_file(file_momenta.nonlocal_current.filename,
			     momenta.nonlocal_current,
			     files.nprop, bar3pt.nrow, first);
	      } // end for n
	    } // end for g
	  } // end for p
	} // end for l
      } // end for f
    } // end for q
  } // end for b
}


void construct_filenames(File_wallformfac_t& files, const WallFormFac_t& bar3pt,
			 int nprop)
{
  files.nprop = nprop;   // the number of configurations

  // Resize to hold how hadronic form-factor blocks are present
  files.bar.resize(bar3pt.bar.size());

  // Loop over hadron contractions
  for(int b=0; b < files.bar.size(); ++b)
  {
    const WallFormFac_bar_t& bar = bar3pt.bar[b];
    File_bar_t& file_bar = files.bar[b];

    // The group of formfactors
    const WallFormFac_formfacs_t& formfacs = bar.formfacs;
    File_formfacs_t& file_formfacs = file_bar.formfacs;

    // Resize to hold how many quark contributions
    file_formfacs.quark.resize(formfacs.quark.size());

    // Loop over the different quark contributions
    for(int q=0; q < file_formfacs.quark.size(); ++q)
    {
      const WallFormFac_quark_t& quark = formfacs.quark[q];
      File_quark_t& file_quark = file_formfacs.quark[q];

      // Resize to hold how many subclasses of formfactors present
      file_quark.formfac.resize(quark.formfac.size());

      // Loop over the different subclasses of formfactors
      for(int f=0; f < file_quark.formfac.size(); ++f)
      {
	const WallFormFac_formfac_t& formfac = quark.formfac[f];
	File_formfac_t& file_formfac = file_quark.formfac[f];

	// Resize to hold how many lorentz terms
	file_formfac.lorentz.resize(formfac.lorentz.size());

	// Loop over the lorentz terms
	for(int l=0; l < file_formfac.lorentz.size(); ++l)
	{
	  const WallFormFac_lorentz_t& lorentz = formfac.lorentz[l];
	  File_lorentz_t& file_lorentz = file_formfac.lorentz[l];

	  // Resize to hold how many projectors
	  file_lorentz.projector.resize(lorentz.projector.size());

	  // Loop over the projectors
	  for(int p=0; p < file_lorentz.projector.size(); ++p)
	  {
	    const WallFormFac_projector_t& projector = lorentz.projector[p];
	    File_projector_t& file_projector = file_lorentz.projector[p];

	    // Resize to hold how many formfactor insertions (gamma matrix values) are writtn
	    file_projector.insertion.resize(projector.insertion.size());

	    // Loop over the gamma matrices
	    for(int g=0; g < file_projector.insertion.size(); ++g)
	    {
	      const WallFormFac_insertion_t& insertion = projector.insertion[g];
	      File_insertion_t& file_insertion = file_projector.insertion[g];

	      // Resize to hold how many momenta values are used
	      file_insertion.momenta.resize(insertion.momenta.size());

	      // Loop over the momenta
	      for(int n = 0; n < file_insertion.momenta.size(); ++n)
	      {
		const WallFormFac_momenta_t& momenta = insertion.momenta[n];
		File_momenta_t& file_momenta = file_insertion.momenta[n];

		ostringstream ff;
		ff << "_" << bar.formfac_type 
		   << "_f" << f
		   << "_" << quark.quark_name 
		   << "_p" << p
		   << "_snk" << lorentz.snk_gamma
		   << "_g" << insertion.gamma_value
		   << "_src" << lorentz.src_gamma
		   << mom_name(momenta.inser_mom); 
      
		file_momenta.local_current.filename = "local_cur3ptfn" + ff.str();
		file_momenta.nonlocal_current.filename = "nonlocal_cur3ptfn" + ff.str();
	      } // end for n
	    } // end for g
	  } // end for p
	} // end for l
      } // end for f
    } // end for q
  } // end for b
}

//
// Main program - loop over files
//
int main(int argc, char **argv)
{
  if (argc < 3)
  {
    cerr << "Usage: " << argv[0] << " <formfac[=all,PION, etc.]>  file1 [file2 file3... fileN]" << endl;
    exit(1);
  }

//  string xml_in_root = "/wallFormFac";

  // The name of the formfac
  string formfac_type = argv[1];

  // Process the first file specially - read alls it params
  int nprop = argc - 2;
  cout << "Reading config " << 1 << " - " << argv[2] << endl;
  BinaryFileReader bin_in(argv[2]);
  
  // Big nested structure that is image of entire file
  WallFormFac_t  wallff;

  // Read data
  clock_t t1 = clock();
  read(bin_in, wallff);
  clock_t t2 = clock();
  cout << "First reading took " << (double)(t2-t1)/(double)(CLOCKS_PER_SEC) << " secs" << endl;

  bin_in.close();

#if 1
  cout << "wallff.bar.size = " 
       << wallff.bar.size() << endl;
  cout << "wallff.bar[0].formfacs.quark.size = " 
       << wallff.bar[0].formfacs.quark.size() << endl;
  cout << "wallff.bar[0].formfacs.quark[0].formfac.size = " 
       << wallff.bar[0].formfacs.quark[0].formfac.size() << endl;
  cout << "wallff.bar[0].formfacs.quark[0].formfac[0].lorentz.size = " 
       << wallff.bar[0].formfacs.quark[0].formfac[0].lorentz.size() << endl;
  cout << "wallff.bar[0].formfacs.quark[0].formfac[0].lorentz[0].projector.size = " 
       << wallff.bar[0].formfacs.quark[0].formfac[0].lorentz[0].projector.size() << endl;
  cout << "wallff.bar[0].formfacs.quark[0].formfac[0].lorentz[0].projector[0].insertion.size = " 
       << wallff.bar[0].formfacs.quark[0].formfac[0].lorentz[0].projector[0].insertion.size() << endl;
  cout << "wallff.bar[0].formfacs.quark[0].formfac[0].lorentz[0].projector[0].insertion[0].momenta.size = " 
       << wallff.bar[0].formfacs.quark[0].formfac[0].lorentz[0].projector[0].insertion[0].momenta.size() << endl;
#endif
  

  // We care about the output version
  switch (wallff.out_version)
  {
  case 3:
  case 4:
    cout << "Processing output version " << wallff.out_version << endl;
    break;

  default:
    cerr << "Unexpected output version " << wallff.out_version << endl;
    throw;
  }


  //
  // Construct file names
  //

  // Structure holding file names
  File_wallformfac_t files;

  // Construct the filenames into the structure  files
  construct_filenames(files, wallff, nprop);

  // Print the first time the files
  t1 = clock();
  cout << "Printing config 1" << endl;
  print_files(formfac_type, files, wallff, true);
  t2 = clock();
  cout << "First printing took " << (double)(t2-t1)/(double)(CLOCKS_PER_SEC) << " secs" << endl;

  // Check filenames
//  cout << "File: " << files.bar[0].seqsrc[4].formfac[1].momenta[0].local_current.filename << endl;

  
  for(int nfile=2; nfile <= nprop; ++nfile)
  {
    // Read data
    /* NOTE: could just read the hadron part */
    cout << "Reading config " << nfile << " - " << argv[nfile+1] << endl;
    bin_in.open(argv[nfile+1]);

    t1 = clock();
    read(bin_in, wallff);
    t2 = clock();
    cout << "Time to read config " << nfile << "  took " << (double)(t2-t1)/(double)(CLOCKS_PER_SEC) << " secs" << endl;


    bin_in.close();

    // Append to the files
    cout << "Printing config " << nfile << endl;
    t1 = clock();
    print_files(formfac_type, files, wallff, false);
    t2 = clock();
    cout << "Time to print config " << nfile << "  took " << (double)(t2-t1)/(double)(CLOCKS_PER_SEC) << " secs" << endl;
  }

  exit(0);
}
