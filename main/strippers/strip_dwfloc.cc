// $Id: strip_dwfloc.cc,v 2.0 2008/12/05 04:44:03 edwards Exp $

// Stripper for the dwf locality test calculation

#include "strippers.h"

using namespace Strippers;
using namespace std;

const int Ns = 4;
const int Nd = 4;

/*
 * Input 
 */

struct Chiral_param_t
{
  Real M5 ;
  int  Ls ;
};

struct Param_t
{
  Array<int>      nrow      ;
  Chiral_param_t  chiral    ;
};

/*
 * Structures for hadron parts
 */

struct Correlator_t
{
  Array<Real> mesprop;
};

struct Measure_t
{
  Correlator_t x_axis ;
  Correlator_t y_axis ;
  Correlator_t z_axis ;
  Correlator_t t_axis ;
};

struct dwfloc_t
{
  Param_t          param          ;
  Measure_t        meas           ;
};



/*
 * Read the parameters
 */

void read(XMLReader& xml, const string& path, Chiral_param_t& chi )
{
  XMLReader top(xml, path);
  read(top, "WilsonMass",chi.M5);
  read(top, "N5", chi.Ls ) ;

}

// Read correlator struct
void read(XMLReader& xml, const string& path, Correlator_t& cor)
{
  XMLReader top(xml, path);
  read(top, "val", cor.mesprop);
}



// Mother of all readers
void read(XMLReader& xml, const string& path, dwfloc_t& dat)
{
  try 
  {
    XMLReader mother(xml, path);

    read(mother, "lattice/size", dat.param.nrow);  // lattice size
    read(mother, "ChiralParam" , dat.param.chiral);
    read(mother, "x_axis"      , dat.meas.x_axis);
    read(mother, "y_axis"      , dat.meas.y_axis);
    read(mother, "z_axis"      , dat.meas.z_axis);
    read(mother, "t_axis"      , dat.meas.t_axis);

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

struct Correlator_files_t
{
  File_t x_axis;
  File_t y_axis;
  File_t z_axis;
  File_t t_axis;
};



struct File_dwfloc_t
{
  int ncor;
  Correlator_files_t  cor;
};



// Write mesons
ostream& operator<<(ostream& s, const Array<Real>& p)
{
  for(int i=0; i < p.size(); ++i)
    s << i << " " << p[i] << "\n";
}

void print_file(const string& filename, const Array<Real>& mesprop, 
		int nprop, 
		const dwfloc_t& dat, int dir,
		bool first)
{
  if (first)
  {
    ofstream file(filename.c_str());

    // Write header line
    file << nprop << " " << dat.param.nrow[dir] << " 0 " 
	 << dat.param.nrow[0] << " 1" << endl;

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


void print_files(const File_dwfloc_t& f, const dwfloc_t& dat, bool first)
{
  
  print_file(f.cor.x_axis.filename,dat.meas.x_axis.mesprop,f.ncor,dat,0,first);
  print_file(f.cor.y_axis.filename,dat.meas.y_axis.mesprop,f.ncor,dat,1,first);
  print_file(f.cor.z_axis.filename,dat.meas.z_axis.mesprop,f.ncor,dat,2,first);
  print_file(f.cor.t_axis.filename,dat.meas.t_axis.mesprop,f.ncor,dat,3,first);
}


void construct_filenames(File_dwfloc_t& files, const dwfloc_t& dat,int nc)
{
  files.ncor = nc;   // the number of configurations
  
  ostringstream lat_ostr ;
  for(int i(0);i<4;i++)
    lat_ostr<<dat.param.nrow[i] ;
  lat_ostr<<"ls"<<dat.param.chiral.Ls ;
  lat_ostr<<"M"<<dat.param.chiral.M5 ;
  
  files.cor.x_axis.filename = "x_axis"+lat_ostr.str() ;
  files.cor.y_axis.filename = "y_axis"+lat_ostr.str() ;
  files.cor.z_axis.filename = "z_axis"+lat_ostr.str() ;
  files.cor.t_axis.filename = "t_axis"+lat_ostr.str() ;
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

  string xml_in_root = "/DWF_locality";


  // Process the first file specially - read alls it params
  cerr << "Open file " << argv[1] << endl;
  XMLReader xml_in(argv[1]);
  cerr << "Done Open file " << argv[1] << endl;
  
  // Big nested structure that is image of entire file
  dwfloc_t  dat;

  // Read data
  read(xml_in, xml_in_root, dat);
  xml_in.close();
  
  cout << " Correlator size = " << dat.meas.x_axis.mesprop.size() << endl;
  
  
  // Structure holding file names
  File_dwfloc_t files;

  // Construct the filenames into the structure  files
  construct_filenames(files, dat, argc-1);

  cout << "x_axis file = ";
  cout << files.cor.x_axis.filename << endl;
  cout << "y_axis file = ";
  cout << files.cor.y_axis.filename << endl;
  cout << "z_axis file = ";
  cout << files.cor.z_axis.filename << endl;
  cout << "t_axis file = ";
  cout << files.cor.t_axis.filename << endl;
  
  // Print the first time the files
  cout << "Printing config 1" << endl;
  print_files(files, dat, true);

  for(int nfile=2; nfile < argc; ++nfile)
  {
    // Read data
    /* NOTE: could just read the hadron part */
    xml_in.open(argv[nfile]);
    read(xml_in, xml_in_root, dat);
    xml_in.close();

    // Append to the files
    cout << "Printing config " << nfile << endl;
    print_files(files, dat, false);
  }

  exit(0);
}
