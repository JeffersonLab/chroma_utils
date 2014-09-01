// $Id: strip_bar3ptfn.cc,v 2.0 2008/12/05 04:44:03 edwards Exp $

#include <iostream>
#include <fstream>

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
  int          j_decay;
  Array<int>   nrow;
};

/*
 * Structures for hadron parts
 */
struct Momenta_t
{
  int            magic;     // magic number for sanity checks
  Array<int>     inser_mom;
  Array<Complex> local_current;
  Array<Complex> nonlocal_current;
};

struct FormFac_t
{
  int              gamma_value;
  Array<Momenta_t> momenta;
};

struct Insertions_t
{
  int               output_version;
  Array<FormFac_t>  formFac;
};

struct Sequential_source_seq_hadron_t
{
  string            seq_src;
  Array<int>        t_srce;
  int               t_sink;
  Array<int>        sink_mom;
  Complex           seq_hadron_0;
  Insertions_t      formFacs;
};

struct Sequential_source_int_t
{
  int               seq_src_value;    // this is now obsolete - convert to string
  Array<int>        t_srce;
  int               t_sink;
  Array<int>        sink_mom;
  Complex           seq_hadron_0;
  Insertions_t      formFacs;
};

struct Sequential_source_t
{
  string            seq_src;
  int               t_source;
  int               t_sink;
  Array<int>        sink_mom;
  int               gamma_insertion;
  Insertions_t      formFacs;

  // Convert old stuff to new format
  Sequential_source_t& operator=(const Sequential_source_seq_hadron_t& s)
  {
    seq_src   = s.seq_src;
    t_source  = s.t_srce[3];
    t_sink    = s.t_sink;
    sink_mom  = s.sink_mom;
    gamma_insertion = 0;
    formFacs  = s.formFacs;
  }

  // Convert old stuff to new format
  Sequential_source_t& operator=(const Sequential_source_int_t& s)
  {
    switch(s.seq_src_value)
    {
    case 0:
      seq_src  = "NUCL_U_UNPOL";
      break;
    case 1:
      seq_src  = "NUCL_D_UNPOL";
      break;
    case 2:
      seq_src  = "NUCL_U_POL";
      break;
    case 3:
      seq_src  = "NUCL_D_POL";
      break;
    case 4:
      seq_src  = "DELTA_U_UNPOL";
      break;
    case 5:
      seq_src  = "DELTA_D_UNPOL";
      break;
    case 6:
      seq_src  = "NUCL_U_UNPOL_NONREL";
      break;
    case 7:
      seq_src  = "NUCL_D_UNPOL_NONREL";
      break;
    case 8:
      seq_src  = "NUCL_U_POL_NONREL";
      break;
    case 9:
      seq_src  = "NUCL_D_POL_NONREL";
      break;
    case 21:
      seq_src  = "NUCL_U_MIXED_NONREL";
      break;
    case 22:
      seq_src  = "NUCL_D_MIXED_NONREL";
      break;
    case 10:
      seq_src  = "PION";
      break;
    default:
      cerr << "Unsupported seq_source_value: " << s.seq_src_value << endl;
      exit(1);
    }

    t_source = s.t_srce[3];
    t_sink = s.t_sink;
    sink_mom = s.sink_mom;
    gamma_insertion = 0;
    formFacs = s.formFacs;
  }
};

struct Wilson_3Pt_fn_measurements_t
{
  int               output_version;   // Unique id for each output version of the structures
  Array<Sequential_source_t> seqsrc;
};

struct Bar3ptfn_t
{
  Output_version_t output_version;
  Param_t          param;
  Wilson_3Pt_fn_measurements_t  bar;
};


/*
 * Similar structure to Bar3ptfn_t that holds filenames
 * NOTE: this cannot be merged into Bar3ptfn_t since that
 * is volatile - upon reading this struct old data is
 * wiped out
 */
/*
 * Structures for hadron parts
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

struct File_formFac_t
{
  Array<File_momenta_t> momenta;
};

struct File_insertions_t
{
  Array<File_formFac_t> formFac;
};

struct File_sequential_source_t
{
  File_insertions_t      formFacs;
};

struct File_wilson_3Pt_fn_measurements_t
{
  Array<File_sequential_source_t> seqsrc;
};

struct File_bar3ptfn_t
{
  int nprop;
  File_wilson_3Pt_fn_measurements_t  bar;
};



// Read formfactors

/*
 * Read the parameters
 */
void read(BinaryReader& bin, Output_version_t& out)
{
  // Not setting context here as it is only the one read.
  read(bin, out.out_version);
}

// params
void read(BinaryReader& bin, Param_t& param)
{
  read(bin, param.mom2_max);
  read(bin, param.j_decay);
  read(bin, param.nrow);
//  cout << "nrow.size() = " << param.nrow.size() << endl;
}


// Read a momenta struct
void read(BinaryReader& bin, Momenta_t& mom)
{
  read(bin, mom.magic);
  if (mom.magic != 20301)
  {
    cerr << "read(Momenta_t): magic number invalid" << endl;
    exit(1);
  }
  read(bin, mom.inser_mom);
  read(bin, mom.local_current);
//  cout << "local_current.size() = " << mom.local_current.size() << endl;
  read(bin, mom.nonlocal_current);
//  cout << "nonlocal_current.size() = " << mom.nonlocal_current.size() << endl;
}

// 
void read(BinaryReader& bin, FormFac_t& mes)
{
  read(bin, mes.gamma_value);
  read(bin, mes.momenta);
//  cout << "momenta.size() = " << mes.momenta.size() << endl;
}

// 
void read(BinaryReader& bin, Insertions_t& mes)
{
  read(bin, mes.output_version);
  read(bin, mes.formFac);
//  cout << "formFac.size() = " << mes.formFac.size() << endl;
}

// 
void read(BinaryReader& bin, Sequential_source_t& src)
{
  read(bin, src.seq_src, 100);
  read(bin, src.t_source);
  read(bin, src.t_sink);
  read(bin, src.sink_mom);
  read(bin, src.gamma_insertion);
  read(bin, src.formFacs);
}

// 
void read(BinaryReader& bin, Sequential_source_seq_hadron_t& src)
{
  read(bin, src.seq_src, 100);
  read(bin, src.t_srce);
  read(bin, src.t_sink);
  read(bin, src.sink_mom);
  read(bin, src.seq_hadron_0);
  read(bin, src.formFacs);
}

// 
void read(BinaryReader& bin, Sequential_source_int_t& src)
{
  read(bin, src.seq_src_value);
  read(bin, src.t_srce);
  read(bin, src.t_sink);
  read(bin, src.sink_mom);
  read(bin, src.seq_hadron_0);
  read(bin, src.formFacs);
}

// Read a hadron measurements
void read(BinaryReader& bin, Wilson_3Pt_fn_measurements_t& had)
{
  read(bin, had.output_version);

  switch (had.output_version)
  {
    case 2:
    {
      Array<Sequential_source_int_t> seqsrc;
      read(bin, seqsrc);
      had.seqsrc.resize(seqsrc.size());
      for(int i=0; i < had.seqsrc.size(); ++i)
	had.seqsrc[i] = seqsrc[i];
    }
    break;

    case 3:
    {
      Array<Sequential_source_seq_hadron_t> seqsrc;
      read(bin, seqsrc);
      had.seqsrc.resize(seqsrc.size());
      for(int i=0; i < had.seqsrc.size(); ++i)
	had.seqsrc[i] = seqsrc[i];

      read(bin, had.seqsrc);
    }
    break;

    case 4:
      read(bin, had.seqsrc);
      break;

    default:
      cerr << "Wilson_3Pt_fn_measurements_t: Unsupported output version " 
	   << had.output_version << endl;
      exit(1);
    }
//  cout << "seqsrc.size() = " << had.seqsrc.size() << endl;
}

// Mother of all readers
void read(BinaryReader& bin, Bar3ptfn_t& bar)
{
  try 
  {
    read(bin, bar.output_version);
    switch (bar.output_version.out_version)
    {
    case 9:
    case 10:
      cerr << "read(Bar3ptfn): oops, need to fix to support old versions" << endl;
      break;

    case 11:
      read(bin, bar.param);
      break;

    default:
      cerr << "read(Bar3ptfn): Unsupported output version " 
	   << bar.output_version.out_version<< endl;
      exit(1);
    }

    read(bin, bar.bar);
  }
  catch( const string& error) { 
    cerr << "Error reading data: " << error << endl;
    exit(1);
  }

  cout << "j_decay=" << bar.param.j_decay;
  cout << " nrow[0]=" << bar.param.nrow[0] 
       << " nrow[1]=" << bar.param.nrow[1] 
       << " nrow[2]=" << bar.param.nrow[2] 
       << " nrow[3]=" << bar.param.nrow[3] 
       << endl;
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


string sink_name(const Array<int>& sink_mom)
{
  ostringstream s;

  // Always print momenta structure
  s << "_pfx" << sink_mom[0] << "_pfy" << sink_mom[1] << "_pfz" << sink_mom[2];

  return s.str();
}


void print_file(const string& filename, const Array<Complex>& current, 
		int nprop, 
		const Param_t& param, 
		bool first)
{
  if (first)
  {
    ofstream file(filename.c_str());

    // Write header line
    // Assuming j_decay != 0
    file << nprop << " " << param.nrow[param.j_decay] << " 1 " 
	 << param.nrow[0] << " 1" << endl;

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


void print_files(const File_bar3ptfn_t& files, const Bar3ptfn_t& bar3pt, bool first)
{
  for(int s=0; s < files.bar.seqsrc.size(); ++s)
    for(int g=0; g < files.bar.seqsrc[s].formFacs.formFac.size(); ++g)
      for(int n = 0; n < files.bar.seqsrc[s].formFacs.formFac[g].momenta.size(); ++n)
      {
	// The local current should always be there, but check anyway
	if (bar3pt.bar.seqsrc[s].formFacs.formFac[g].momenta[n].local_current.size() > 0)
	  print_file(files.bar.seqsrc[s].formFacs.formFac[g].momenta[n].local_current.filename,
		     bar3pt.bar.seqsrc[s].formFacs.formFac[g].momenta[n].local_current,
		     files.nprop, bar3pt.param, first);

	// The non-local current may not always be present, so check
	if (bar3pt.bar.seqsrc[s].formFacs.formFac[g].momenta[n].nonlocal_current.size() > 0)
	  print_file(files.bar.seqsrc[s].formFacs.formFac[g].momenta[n].nonlocal_current.filename,
		     bar3pt.bar.seqsrc[s].formFacs.formFac[g].momenta[n].nonlocal_current,
		     files.nprop, bar3pt.param, first);
      }
}


void construct_filenames(File_bar3ptfn_t& files, const Bar3ptfn_t& bar3pt,
			 int nprop)
{
  files.nprop = nprop;   // the number of configurations

  // Resize to hold how many seq. sources are involved
  files.bar.seqsrc.resize(bar3pt.bar.seqsrc.size());

  // Loop over seq. srcs
  for(int s=0; s < files.bar.seqsrc.size(); ++s)
  {
    string seq_src = bar3pt.bar.seqsrc[s].seq_src;
    int gamma_insertion = bar3pt.bar.seqsrc[s].gamma_insertion;    // Gamma insertion

    // Resize to hold how many formfactor insertions (gamma matrix values) are writtn
    files.bar.seqsrc[s].formFacs.formFac.resize(bar3pt.bar.seqsrc[s].formFacs.formFac.size());

    // Loop over the gamma matrices
    for(int g=0; g < files.bar.seqsrc[s].formFacs.formFac.size(); ++g)
    {
      int gamma_value = bar3pt.bar.seqsrc[s].formFacs.formFac[g].gamma_value;

      // Resize to hold how many momenta values are used
      files.bar.seqsrc[s].formFacs.formFac[g].momenta.resize(
	bar3pt.bar.seqsrc[s].formFacs.formFac[g].momenta.size());

      // Loop over the momenta
      for(int n = 0; n < files.bar.seqsrc[s].formFacs.formFac[g].momenta.size(); ++n)
      {
	ostringstream f;
	f << "_" << seq_src << "_i" << gamma_insertion << "_g" << gamma_value
	  << mom_name(bar3pt.bar.seqsrc[s].formFacs.formFac[g].momenta[n].inser_mom)
	  << sink_name(bar3pt.bar.seqsrc[s].sink_mom);

	files.bar.seqsrc[s].formFacs.formFac[g].momenta[n].local_current.filename =
	  "local_cur3ptfn" + f.str();

	files.bar.seqsrc[s].formFacs.formFac[g].momenta[n].nonlocal_current.filename =
	  "nonlocal_cur3ptfn" + f.str();
      } // end for n
    } // end for g
  } // end for s
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

  // Process the first file specially - read alls it params
  clock_t t1 = clock();
  BinaryFileReader bin_in(argv[1]);
  clock_t t2 = clock();
  cout << "First read took " << (double)(t2-t1)/(double)(CLOCKS_PER_SEC) << " secs" << endl;
  
  // Big nested structure that is image of entire file
  Bar3ptfn_t  bar3pt;

  // Read data
  t1 = clock();
  read(bin_in, bar3pt);
  t2 = clock();
  cout << "First processing took " << (double)(t2-t1)/(double)(CLOCKS_PER_SEC) << " secs" << endl;

  bin_in.close();

  cout << "bar3pt.bar.seqsrc.size = " 
       << bar3pt.bar.seqsrc.size() << endl;
  cout << "bar3pt.bar.seqsrc.formFacs.formFac.size = " 
       << bar3pt.bar.seqsrc[0].formFacs.formFac.size() << endl;
  cout << "bar3pt.bar.seqsrc.formFacs.formFac[0].momenta.size = " 
       << bar3pt.bar.seqsrc[0].formFacs.formFac[0].momenta.size() << endl;
  cout << "bar3pt.bar.seqsrc.formFacs.formFac[0].momenta.local_current.size = " 
       << bar3pt.bar.seqsrc[0].formFacs.formFac[0].momenta[0].local_current.size() << endl;
  

  // We care about the output version
  switch (bar3pt.output_version.out_version)
  {
  case 9:
  case 10:
  case 11:
    cout << "Processing output version " << bar3pt.output_version.out_version << endl;
    break;

  default:
    cerr << "Unexpected output version " << bar3pt.output_version.out_version << endl;
    exit(1);
  }


  //
  // Construct file names
  //

  // Structure holding file names
  File_bar3ptfn_t files;

  // Construct the filenames into the structure  files
  construct_filenames(files, bar3pt, argc-1);

  // Print the first time the files
  t1 = clock();
  cout << "Printing config 1" << endl;
  print_files(files, bar3pt, true);
  t2 = clock();
  cout << "First printing took " << (double)(t2-t1)/(double)(CLOCKS_PER_SEC) << " secs" << endl;

  // Check filenames
//  cout << "File: " << files.bar.seqsrc[4].formFacs.formFac[1].momenta[0].local_current.filename << endl;

  
  for(int nfile=2; nfile < argc; ++nfile)
  {
    // Read data
    /* NOTE: could just read the hadron part */
    cout << "Reading config " << nfile << " - " << argv[nfile] << endl;
    t1 = clock();
    bin_in.open(argv[nfile]);
    t2 = clock();
    cout << "Time to read config " << nfile << "  took " << (double)(t2-t1)/(double)(CLOCKS_PER_SEC) << " secs" << endl;

    t1 = clock();
    read(bin_in, bar3pt);
    t2 = clock();
    cout << "Time to process config " << nfile << "  took " << (double)(t2-t1)/(double)(CLOCKS_PER_SEC) << " secs" << endl;
    bin_in.close();

    // Append to the files
    cout << "Printing config " << nfile << endl;
    t1 = clock();
    print_files(files, bar3pt, false);
    t2 = clock();
    cout << "Time to print config " << nfile << "  took " << (double)(t2-t1)/(double)(CLOCKS_PER_SEC) << " secs" << endl;
  }

  exit(0);
}
