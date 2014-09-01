// $Id: strip_buildingblocks.cc,v 2.0 2008/12/05 04:44:03 edwards Exp $

#include <iostream>
#include <fstream>

#include <sys/time.h>   // for timings
#include <assert.h>

#include "strippers.h"
#include "io/adat_byteorder.h"
#include "formfac/formfac.h"

using namespace Strippers;
using namespace std;
using namespace ENSEM;
using namespace FF;

/*
 * Similar structure to BuildingBlocks_t that holds filenames
 * NOTE: this cannot be merged into BuildingBlocks_t since that
 * is volatile - upon reading this struct old data is
 * wiped out
 */
/*
 * Structures for hadron parts
 */
struct File_BuildingBlocks_t
{
  struct File_links_t
  {
    struct File_gamma_t
    {
      struct File_t
      {
	std::string    filename;
      };

      Array<File_t>          momenta; 
    };

    Array<File_gamma_t>      gamma;
  };
  
  int                        nprop;
  Array<File_links_t>        links;
};



// Write formfactors
ostream& operator<<(ostream& s, const Array<ComplexF>& p)
{
  for(int i=0; i < p.size(); ++i)
    s << i << " " << real(p[i]) << " " << imag(p[i]) << "\n";
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



string link_name(const Array<int>& links)
{
  ostringstream s;

  s << "_l" << links.size();

  for(int i=0; i < links.size(); ++i)
  {
//    s << "_" << ((links[i] < 0) ? '-' : '+') << links[i];
    s << "_" << links[i];
  }

  return s.str();
}


string mom_name(const Array<int>& sink_mom)
{
  ostringstream s;

  // Always print momenta structure
  s << "_qx" << sink_mom[0] << "_qy" << sink_mom[1] << "_qz" << sink_mom[2];

  return s.str();
}


void print_file(const string& filename, const Array<ComplexF>& current, 
		int nprop, 
		const BuildingBlocks_t::Param_t& param, 
		bool first)
{
  if (first)
  {
    ofstream file(filename.c_str());

    // Write header line
// FIX THIS - RATS - need to GET j_decay BACK IN HERE!!
    file << nprop << " " << param.latt_size[Nd-1] << " 1 " 
	 << param.latt_size[0] << " 1" << endl;

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


void print_files(const File_BuildingBlocks_t& files, const BuildingBlocks_t& bar3pt, bool first)
{
  for(int s=0; s < files.links.size(); ++s)
    for(int g=0; g < files.links[s].gamma.size(); ++g)
      for(int m=0; m < files.links[s].gamma[g].momenta.size(); ++m)
      {
	// The local current should always be there, but check anyway
	print_file(files.links[s].gamma[g].momenta[m].filename,
		   bar3pt.links[s].gamma[g].momenta[m].current,
		   files.nprop, bar3pt.param, first);
      }
}


void construct_filenames(File_BuildingBlocks_t& files, const BuildingBlocks_t& bar3pt,
			 int nprop, const string& seq_src)
{
  files.nprop = nprop;   // the number of configurations


  // Resize to hold how many seq. sources are involved
  files.links.resize(bar3pt.links.size());

  // Loop over links
  for(int l=0; l < files.links.size(); ++l)
  {
    // Resize to hold how many formfactor insertions (gamma matrix values) are writtn
    files.links[l].gamma.resize(bar3pt.links[l].gamma.size());

    string ss = link_name(bar3pt.links[l].link_value);

    // Loop over the gamma matrices
    for(int g=0; g < files.links[l].gamma.size(); ++g)
    {
      // Resize to hold number of momenta
      files.links[l].gamma[g].momenta.resize(bar3pt.links[l].gamma[g].momenta.size());

      int gamma_value = bar3pt.links[l].gamma[g].gamma_value;
      
      // Loop over momenta
      for(int m=0; m < files.links[l].gamma[g].momenta.size(); ++m)
      {
	Array<int> inser_mom = bar3pt.links[l].gamma[g].momenta[m].inser_mom;

	ostringstream f;
	f << ss << "_g" << gamma_value << mom_name(inser_mom); 

	files.links[l].gamma[g].momenta[m].filename = seq_src + f.str();
      } // end for m
    } // end for g
  } // end for s
}



//
// Main program - loop over files
//
int main(int argc, char **argv)
{
  if (argc <= 1)
  {
    cerr << "Usage: " << argv[0] << " PREFIX file1 [file2 file3... fileN]" << endl;
    exit(1);
  }

  // The prefix of the stripped files
  // THIS SHOULD REALLY BE THE SEQSOURCE TYPE, BUT BB DOESN'T RECORD IT
  const string prefix = argv[1];

  // Check if the prefix really is a file. That is a mistake
  {
    FILE* fp = fopen(prefix.c_str(), "rb");
    if (fp != 0)
    {
      fclose(fp);
      cerr << "Usage: " << argv[0] << " PREFIX file1 [file2 file3... fileN]" << endl;
      exit(1);
    }
  }

  // Process the first file specially - read alls it params
  int nprop = argc - 2;
  cout << "Reading config " << 1 << " - " << argv[2] << endl;
  
  // Big nested structure that is image of entire file
  BuildingBlocks_t  bar3pt;

  // Read data
  clock_t t1 = clock();
  read(argv[2], bar3pt);
  clock_t t2 = clock();
  cout << "First processing took " << (double)(t2-t1)/(double)(CLOCKS_PER_SEC) << " secs" << endl;

  cout << "bar3pt.links.size = " 
       << bar3pt.links.size() << endl;
  cout << "bar3pt.links[0].gamma.size = " 
       << bar3pt.links[0].gamma.size() << endl;
  cout << "bar3pt.links.gamma[0].momenta.size = " 
       << bar3pt.links[0].gamma[0].momenta.size() << endl;
  cout << "bar3pt.links.gamma[0].momenta.current.size = " 
       << bar3pt.links[0].gamma[0].momenta[0].current.size() << endl;
  
  //
  // Construct file names
  //
  // Structure holding file names
  File_BuildingBlocks_t files;

  // Construct the filenames into the structure  files
  construct_filenames(files, bar3pt, nprop, prefix);

  // Print the first time the files
  t1 = clock();
  cout << "Printing config 1" << endl;
  print_files(files, bar3pt, true);
  t2 = clock();
  cout << "First printing took " << (double)(t2-t1)/(double)(CLOCKS_PER_SEC) << " secs" << endl;

  
  for(int nfile=2; nfile <= nprop; ++nfile)
  {
    // Read data
    /* NOTE: could just read the hadron part */
    cout << "Reading config " << nfile << " - " << argv[nfile+1] << endl;
    t1 = clock();
    read(argv[nfile+1], bar3pt);
    t2 = clock();
    cout << "Time to process config " << nfile << "  took " << (double)(t2-t1)/(double)(CLOCKS_PER_SEC) << " secs" << endl;

    // Append to the files
    cout << "Printing config " << nfile << endl;
    t1 = clock();
    print_files(files, bar3pt, false);
    t2 = clock();
    cout << "Time to print config " << nfile << "  took " << (double)(t2-t1)/(double)(CLOCKS_PER_SEC) << " secs" << endl;
  }

  exit(0);
}
