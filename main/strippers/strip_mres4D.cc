// $Id: strip_mres4D.cc,v 2.0 2008/12/05 04:44:04 edwards Exp $

// Stripper for the Mres and Za calculation

#include "strippers.h"

using namespace Strippers;
using namespace std;


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

  string xml_in_root;
  string xml_in_root_a = "/mres";
  string xml_in_root_b = "/t_mres4D";

  int n_corr=argc-1; 
  cout << "Trying to read " << n_corr << " files " << endl;

  int n_t;

  Array< Array<double> > delta_corr(n_corr);
  Array< Array<double> > delta_sq_corr(n_corr);
  Array< Array<double> > pseudo_corr(n_corr);


  for(int i=0; i < n_corr; i++) { 
    string filename_in(argv[i+1]);

    cout << "File is: " << filename_in << endl;
    try { 
      XMLReader xml_in(filename_in);
      
      if( xml_in.count(xml_in_root_a)==1 )
	xml_in_root = xml_in_root_a;
      else if( xml_in.count(xml_in_root_b)==1 )
	xml_in_root = xml_in_root_b;
      else
	throw string("unknown root tag");

      read( xml_in, xml_in_root+"/DeltaProp_correlator/delta_prop_corr", delta_corr[i]);
      read( xml_in, xml_in_root+"/DeltaSqProp_correlator/delta_sq_prop_corr", delta_sq_corr[i]);


      if( xml_in.count(xml_in_root+"/PsuedoPseudo_correlator/psudo_prop_corr")==1 ) {
	// Backward compatibility mode -- pre typo days
      read( xml_in, xml_in_root+"/PsuedoPseudo_correlator/psudo_prop_corr", pseudo_corr[i]);
      }
      else {
	 read( xml_in, xml_in_root+"/PsuedoPseudo_correlator/pseudo_prop_corr", pseudo_corr[i]);
      }
    }
    catch(const string& e) { 
      cerr << "Caught Exception: " << e << endl;
      abort();
    }

    if(i==0) { 
      n_t = delta_corr[0].size();
    }
    
    if( delta_corr[i].size() != n_t ) { 
      cerr << "Delta Corr["<<i<<"] has wrong number of elements " << endl;
      abort();
    }

    if( pseudo_corr[i].size() != n_t ) { 
      cerr << "Pseudo Corr["<<i<<"] has wrong number of elements " << endl;
      abort();
    }
  }

  string delta_file_out("./delta4d_corr.dat");
  string deltasq_file_out("./deltaSq4d_corr.dat");
  string pseudo_file_out("./pseudo4d_corr.dat");
  
  ofstream delta_file(delta_file_out.c_str());
  ofstream deltasq_file(deltasq_file_out.c_str());
  ofstream pseudo_file(pseudo_file_out.c_str());

  delta_file << n_corr << " " << n_t << " 0 1 1" << endl;
  deltasq_file << n_corr << " " << n_t << " 0 1 1" << endl;
  pseudo_file << n_corr << " " << n_t << " 0 1 1" << endl;

  for(int i=0; i < n_corr ; i++) { 
    for(int t=0; t < n_t; t++) { 
      delta_file << t << " " << (delta_corr[i])[t] << endl;
      deltasq_file << t << " " << (delta_sq_corr[i])[t] << endl;
      pseudo_file << t << " " << (pseudo_corr[i])[t] << endl;

    }
  }

  delta_file.close();
  deltasq_file.close();
  pseudo_file.close();

  ofstream sum_delta_file("./sum_delta4d.dat");
  ofstream sum_deltasq_file("./sum_deltasq4d.dat");
  ofstream sum_pseudo_file("./sum_pseudo4d.dat");

  sum_delta_file << n_corr << " 1 0 0 1 " << endl;
  sum_deltasq_file << n_corr << " 1 0 0 1 " << endl;
  sum_pseudo_file << n_corr << " 1 0 0 1 " << endl;

  Array<double> sum_delta_c(n_corr);
  Array<double> sum_deltasq_c(n_corr);
  double ensemble_ave_pseudo=0;

  for(int i=0; i < n_corr; i++) { 
    double sum_delta=0;
    double sum_deltasq=0;
    double sum_pseudo=0;

    // Sum over time
    for(int t=0; t < n_t; t++) { 
      sum_delta += (delta_corr[i])[t];
      sum_deltasq += (delta_sq_corr[i])[t];
      sum_pseudo += (pseudo_corr[i])[t];
    }
    sum_delta_file << "0 " << sum_delta << endl;
    sum_deltasq_file << "0 " << sum_deltasq << endl;
    sum_pseudo_file << "0 " << sum_pseudo << endl;
    sum_delta_c[i] = sum_delta;
    sum_deltasq_c[i] = sum_deltasq;
    ensemble_ave_pseudo += sum_pseudo;
  }
  sum_delta_file.close();
  sum_deltasq_file.close();
  sum_pseudo_file.close();

  ensemble_ave_pseudo /= (double)n_corr;
  
  ofstream mres_per_conf_file("./mres4d_per_conf.dat");
  ofstream delsq_per_conf_file("./mres4d_DeltaSq_per_conf.dat");

  for(int i=0; i < n_corr; i++) { 
    mres_per_conf_file << i+1 << " " << sum_delta_c[i] / ensemble_ave_pseudo << endl;
    delsq_per_conf_file << i+1 << " " << sum_deltasq_c[i] / ensemble_ave_pseudo << endl;
  }
  mres_per_conf_file.close();
  delsq_per_conf_file.close();


  double mres=0;
  double mressq=0;
  for(int i=0;i < n_corr; i++) { 
    mres += sum_delta_c[i] / ensemble_ave_pseudo;
    mressq += sum_deltasq_c[i] /ensemble_ave_pseudo;
  }
  mres /= (double)n_corr;
  mressq /= (double)n_corr;

  std::cout << "M_res (guess without errors) = " << mres << endl;
  std::cout << "Delta_Sq (guess without errors) = " << mressq << endl;

  exit(0);
}
