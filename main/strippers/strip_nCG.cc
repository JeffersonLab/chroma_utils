// $Id: strip_nCG.cc,v 2.0 2008/12/05 04:44:04 edwards Exp $

// Stripper for the Mres and Za calculation

#include "strippers.h"

using namespace Strippers;
using namespace std;

using XMLArray::Array;

const int Ns = 4;
const int Nd = 4;

/*
 * Input 
 */

// Main program - loop over files
//
int main(int argc, char **argv)
{
  if (argc == 1)
  {
    cerr << "Usage: " << argv[0] << " file1 [file2 file3... fileN]" << endl;
    exit(1);
  }

  string xml_in_root="//Relaxation_Iterations";

  int n_corr=argc-1; 
  cout << "Trying to read " << n_corr << " files " << endl;

  int n_t;

  Array< double > ncg_had(n_corr);


  for(int i=0; i < n_corr; i++) { 
    string filename_in(argv[i+1]);

    cout << "File is: " << filename_in << endl;
    try { 
      XMLReader xml_in(filename_in);
      
      read( xml_in, xml_in_root+"/ncg_had", ncg_had[i]);
  
    }
    catch(const string& e) { 
      cerr << "Caught Exception: " << e << endl;
      abort();
    }

  }

  string ncg_file_out("./ncg_had.dat");
  ofstream ncg_file(ncg_file_out.c_str());

  ncg_file << n_corr << " " << 1 << " 0 0 1" << endl;

  for(int i=0; i < n_corr ; i++) { 
      ncg_file << "0" << " " << ncg_had[i] << endl;
  }

  ncg_file.close();
  exit(0);
}
