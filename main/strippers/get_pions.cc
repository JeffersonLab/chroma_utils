
#include <iostream>
#include <sstream>

using namespace std;

#include "ensem/ensem.h"
#include "xml_array.h"
#include "io/adat_xmlio.h"
#include "io/adat_io.h"

using XMLArray::Array;
using namespace ENSEM;

struct Pion { 
  double mass;
  Array< Array<double > > corr_samples; 
};

void read_dwf_pion(Pion& pion, const Array<string>& filenames);
void read_ovlap_pions(Array<Pion>& pion, const Array<string>& filenames);

int main(int argc, char *argv[]) { 


  // Count the number of command line args - the command
  if( argc == 1 ) { 
    cerr << "Usage: get_pions <xml_file1> <xml_file2> ... " << endl;
    exit(1);
  }

  int n_input_files=argc-1;
  Array<string> input_filenames(n_input_files);
  cout << "No of Input Files=" << n_input_files <<endl;
  for(int i=0; i < n_input_files; i++) { 
    input_filenames[i] = string(argv[i+1]);
  }
  cout << "Input filenames:" << endl ;
  for(int i=0; i < input_filenames.size(); i++) { 
    cout << "   " << input_filenames[i] << endl;
  }

  // Read the data cubes from the XML files (cleaning)
  Pion dwf_pion;
  Array<Pion> ovlap_pions;

  read_dwf_pion(dwf_pion, input_filenames);
  read_ovlap_pions(ovlap_pions, input_filenames);

  // Transfer to data cubes 
  int n_samples = dwf_pion.corr_samples.size();
  if( n_samples == 0 ) { 
    cerr << "No pion samples are available" << endl;
    exit(1);
  }
  int sample_size = dwf_pion.corr_samples[0].size();
  if( sample_size == 0 ) { 
    cerr << "The DWF sample size appears to be zero" << endl;
  }

  EnsemVectorReal dwf_pions_ens;
  EnsemVectorReal pions_ens;
  EnsemVectorReal ratio_ens;

  dwf_pions_ens.resize(n_samples);
  dwf_pions_ens.resizeObs(sample_size);

  pions_ens.resize(n_samples);
  pions_ens.resizeObs(sample_size);

  ratio_ens.resize(n_samples);
  ratio_ens.resizeObs(sample_size);


  // Get the DWF pions
  for(int i=0; i < n_samples; i++) { 
    for(int j=0; j < sample_size; j++) { 
      dwf_pions_ens.elem(i).elem(j) = dwf_pion.corr_samples[i][j];
    }
  }

  // Dump out the ensemble for use in Analysis in "stripped format"
  TextWriter dwf_pion_writer("./dwf_pion_ensemble.ens");
  dwf_pion_writer << dwf_pions_ens ;

  // Do each mass
  for(int m=0; m < ovlap_pions.size(); m++) {

    // Fill out the pion ensemble
    for(int i=0; i < n_samples; i++) {
      for(int j=0; j < sample_size; j++) { 
	pions_ens.elem(i).elem(j) = ovlap_pions[m].corr_samples[i][j];
      }
    }

    // Dump it
    ostringstream pion_ensemble_filename;
    pion_ensemble_filename << "./ovlap_pion_m"<<ovlap_pions[m].mass<<".ens";
    TextWriter ovlap_pion_writer(pion_ensemble_filename.str());
    ovlap_pion_writer << pions_ens;

    // Now do the ratio
    ratio_ens = pions_ens / dwf_pions_ens;
    ostringstream ratio_filename;
    ratio_filename << "./ovlap_pion_dwf_ratio_m"<<ovlap_pions[m].mass <<".out";
    ofstream ratio_file(ratio_filename.str().c_str());
    calc(ratio_file, ratio_ens);

  }

	       

}

// Grok the DWF pion from the input XML files
void read_dwf_pion(Pion& dwf_pion, const Array<string>& input_filenames) {
  int n_samples=input_filenames.size();
  int nt;
  double dwf_mass;
  try { 
    Array<int> nrow;
    XMLReader first_file(input_filenames[0]);
    read(first_file, "/spectrum_w/Input/Param/nrow", nrow);
    cout << "No of tslices=" << nrow[3] << endl;
    read(first_file, "/spectrum_w//Wilson_hadron_measurements//FermionAction[string(FermAct)='DWF']/Mass", dwf_mass);
    cout << "DWF Fermion mass is " << dwf_mass << endl;
    nt = nrow[3];

    for(int i=1; i < n_samples; i++) {
      XMLReader cur_file(input_filenames[i]);
      read(cur_file, "/spectrum_w/Input/Param/nrow", nrow);
      if( nt != nrow[3] ) { 
	cerr << "Nt mismatch. Nt=" << nt << " read " << nrow[3] << " in file " << input_filenames[i] << endl;
	exit(1);
      }
    }
  }
  catch(const string& e) { 
    cerr << "Caught exception while reading XML " << endl;
    exit(1);
  }

  // At this point we have sanity checked that each file has the
  // same number of timeslices. Now resize the ensem vector appropriately
  // n_samples * n_t
  dwf_pion.corr_samples.resize(n_samples);

  // Now fill out the samples
  try { 
    for(int i=0; i < n_samples; i++) { 
      XMLReader curr_file(input_filenames[i]);
      read(curr_file, "/spectrum_w//Wilson_hadron_measurements//FermionAction[string(FermAct)='DWF']/../../Forward_prop_correlator/forward_prop_corr", dwf_pion.corr_samples[i]);
      cout << "Read Correlator: " << endl;
      for(int t=0; t < nt; t++) { 
	cout << dwf_pion.corr_samples[i][t] << " " ;
      }
      cout << endl;
      dwf_pion.mass=dwf_mass;
    }
  }
  catch( const string& e ) { 
    cerr << "Caught exception while reading XML " << endl;
    exit(1);
  }
}

void read_ovlap_pions(Array<Pion>& pion, const Array<string>& filenames)
{
  int n_samples = filenames.size();
  int n_masses;
  int nt;

  try{
    {
      // Open the first file
      XMLReader first_file(filenames[0]);
      
      // Get the number of masses
      n_masses=first_file.count("//Wilson_hadron_measurements/elem//FermionAction[string(FermAct)='OVERLAP_CONTINUED_FRACTION_5D']");
      cout << "No of Overlap Masses = " << n_masses << endl;
      
      pion.resize(n_masses);
      // Now read the masses
      for(int m=0; m < n_masses; m++) {
	ostringstream query;
	query << "(//Wilson_hadron_measurements/elem//FermionAction[string(FermAct)='OVERLAP_CONTINUED_FRACTION_5D']/Mass)[" << m+1 << "]";
	read(first_file,query.str() , pion[m].mass);
	cout << "Mass[" << m << "] = " << pion[m].mass << endl;
      }
      
      // Resize the number of samples for each mass
      for(int m=0; m < n_masses; m++) {
	pion[m].corr_samples.resize(n_samples);
      }
    }

    for(int i=0; i < n_samples; i++) { 
      XMLReader curr_xml(filenames[i]);

      for(int m=0; m < n_masses; m++) {
	ostringstream query;
	query << "(//Wilson_hadron_measurements/elem//FermionAction[string(FermAct)='OVERLAP_CONTINUED_FRACTION_5D']/Mass)[" << m+1 << "]/../../../Forward_prop_correlator/forward_prop_corr";
	read(curr_xml, query.str(), pion[m].corr_samples[i]);
	cout << "File: " << filenames[i] << " m=" << m << " Mass = " << pion[m].mass << endl;
	for(int t=0; t < pion[m].corr_samples[i].size(); t++) { 
	  cout << pion[m].corr_samples[i][t] << " ";
	}
	cout << endl;
      }
    }
     

  }
  catch(const string& e) {
    cerr << "Caught exception reading XML: " << e << endl;
    exit(1);
  }
}
