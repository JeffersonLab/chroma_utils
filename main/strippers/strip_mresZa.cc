// $Id: strip_mresZa.cc,v 2.0 2008/12/05 04:44:04 edwards Exp $

// Stripper for the Mres and Za calculation

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

struct Chiral_param_t
{
  Real M5 ;
  int  Ls ;
};

struct Param_t
{
  int             version   ;
  Array<int>      nrow      ;
  Chiral_param_t  chiral    ;
  string          FermTypeP ;
  Real            Mass      ;
  int             t_dir     ;
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
  Correlator_t conAx ;
  Correlator_t locAx ;
  Correlator_t mp    ;
  Correlator_t ps    ;
};

struct Mres_t
{
  Output_version_t output_version ;
  Param_t          param          ;
  Measure_t        meas           ;
};



/*
 * Read the parameters
 */
void read(XMLReader& xml, const string& path, Output_version_t& out)
{
  read(xml, path + "/out_version", out.out_version) ;
}

void read(XMLReader& xml, const string& path, Chiral_param_t& chi )
{
  XMLReader top(xml, path);
  read(top, "OverMass",chi.M5);
  read(top, "N5", chi.Ls ) ;

}
// params
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader top(xml, path);

  read(top, "version", param.version); 
  switch (param.version){
  case 4:
    read(top, "ChiralParam",param.chiral);
    read(top, "Mass", param.Mass);
    break ;
  case 5:
    read(top, "FermionAction/ChiralParam",param.chiral);
    read(top, "FermionAction/Mass", param.Mass);
    break ;
  case 6:
  case 7:
  case 8:
  case 9:
    read(top, "FermionAction",param.chiral);
    read(top, "FermionAction/Mass", param.Mass);
    break ;

  default:
    cerr << "read(mresZa): param  version " 
	 << param.version<< endl;
    throw;
  }
  
   if(top.count("nrow")!=0)
     read(top, "nrow", param.nrow);  // lattice size
 
  //read(top, "FermTypeP", param.FermTypeP);
  param.FermTypeP = "WILSON" ;

  //cout<<"I am HERE"<<endl ;
}

// Read correlator struct
void read(XMLReader& xml, const string& path, Correlator_t& cor)
{
  XMLReader top(xml, path);
  read(top, "mesprop", cor.mesprop);
}



// Read Mres measurements
void read(XMLReader& xml, const string& path, Measure_t& meas)
{
//  cout << "Wilson_had = " << path << endl;

  XMLReader top(xml, path);


  read(top, "DWF_ConservedAxial"  , meas.conAx);
  read(top, "DWF_LocalAxial"      , meas.locAx);

  read(top, "DWF_MidPoint_Pseudo" , meas.mp   );
  read(top, "DWF_Psuedo_Pseudo"   , meas.ps   );

}

// Mother of all readers
void read(XMLReader& xml, const string& path, Mres_t& mres)
{
  try 
  {
    XMLReader mother(xml, path);

    read(mother, "Output_version", mres.output_version);

    if(mother.count("ProgramInfo/Setgeom/latt_size")!=0)
      read(mother, "ProgramInfo/Setgeom/latt_size", mres.param.nrow);
    else
      cerr<<"ProgramInfo/Setgeom/latt_size not fount!\n" ;

    switch (mres.output_version.out_version)
    {
    case 1:
      
      if (mother.count("Input/propagator/Param") != 0)
        read(mother, "Input/propagator/Param", mres.param);
      else if (mother.count("Input/propagator/param") != 0)
        read(mother, "Input/propagator/param", mres.param);
      else
        read(mother, "Input/Param", mres.param);
      

      //read(mother, "Input/Param", mres.param);
      //read(mother, "Input/propagator/Param", mres.param);
      //read(mother, "Input/propagator/param", mres.param);
      
      read(mother, "DWF_QuarkProp4/time_direction/t_dir", mres.param.t_dir);
      read(mother, "DWF_QuarkProp4",mres.meas);
      break;


    default:
      cerr << "read(mresZa): Unsupported output version " 
	   << mres.output_version.out_version<< endl;
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

struct Correlator_files_t
{
  File_t conAx;
  File_t locAx;
  File_t mp   ;
  File_t ps   ;
};



struct File_Mres_t
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
		const Mres_t& mres, 
		bool first)
{
  if (first)
  {
    ofstream file(filename.c_str());

    // Write header line
    file << nprop << " " << mres.param.nrow[mres.param.t_dir] << " 0 " 
	 << mres.param.nrow[0] << " 1" << endl;

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


void print_files(const File_Mres_t& f, const Mres_t& mres, bool first)
{
  
  print_file(f.cor.conAx.filename,mres.meas.conAx.mesprop,f.ncor,mres,first);
  print_file(f.cor.locAx.filename,mres.meas.locAx.mesprop,f.ncor,mres,first);
  print_file(f.cor.mp   .filename,mres.meas.mp   .mesprop,f.ncor,mres,first);
  print_file(f.cor.ps   .filename,mres.meas.ps   .mesprop,f.ncor,mres,first);
}


void construct_filenames(File_Mres_t& files, const Mres_t& mres,int nc)
{
  files.ncor = nc;   // the number of configurations
  
  ostringstream lat_ostr ;
  for(int i(0);i<4;i++)
    lat_ostr<<mres.param.nrow[i] ;
  lat_ostr<<"ls"<<mres.param.chiral.Ls ;
  lat_ostr<<"M"<<mres.param.chiral.M5 ;
  lat_ostr<<"m"<<mres.param.Mass ;
  
  files.cor.conAx.filename = "conAx"+lat_ostr.str() ;
  files.cor.locAx.filename = "locAx"+lat_ostr.str() ;
  files.cor.mp   .filename = "midps"+lat_ostr.str() ;
  files.cor.ps   .filename = "pseud"+lat_ostr.str() ;
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

  string xml_in_root = "/propagator";


  // Process the first file specially - read alls it params
  cerr << "Open file " << argv[1] << endl;
  XMLReader xml_in(argv[1]);
  cerr << "Done Open file " << argv[1] << endl;
  
  // Big nested structure that is image of entire file
  Mres_t  mres;

  // Read data
  read(xml_in, xml_in_root, mres);
  xml_in.close();
  
  cout << " Correlator size = " << mres.meas.conAx.mesprop.size() << endl;
  
  // We care about the output version
  switch (mres.output_version.out_version)
  {
  case 1:
    cout<<"Processing output version "<<mres.output_version.out_version<<endl;
    break;

  default:
    cerr<<"Unexpected output version "<<mres.output_version.out_version<<endl;
    throw;
  }

  // Check fermion type
  if (mres.param.FermTypeP != "WILSON")
  {
    cerr << "Expecting WILSON fermion type";
    throw;
  }


  // Structure holding file names
  File_Mres_t files;

  // Construct the filenames into the structure  files
  construct_filenames(files, mres, argc-1);

  cout << "Conserved axial current file = ";
  cout << files.cor.conAx.filename << endl;
  cout << "Loxal axial current file = ";
  cout << files.cor.locAx.filename << endl;
  cout << "Midpoint - Pseudoscalar file = ";
  cout << files.cor.mp.filename << endl;
  cout << "Pseudoscalar - Pseudoscalar file = ";
  cout << files.cor.ps.filename << endl;

  // Print the first time the files
  cout << "Printing config 1" << endl;
  print_files(files, mres, true);

  for(int nfile=2; nfile < argc; ++nfile)
  {
    // Read data
    /* NOTE: could just read the hadron part */
    xml_in.open(argv[nfile]);
    read(xml_in, xml_in_root, mres);
    xml_in.close();

    // Append to the files
    cout << "Printing config " << nfile << endl;
    print_files(files, mres, false);
  }

  exit(0);
}
