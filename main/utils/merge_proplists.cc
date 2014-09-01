// $Id: merge_proplists.cc,v 2.0 2008/12/05 04:44:06 edwards Exp $
// $Log: merge_proplists.cc,v $
// Revision 2.0  2008/12/05 04:44:06  edwards
// Changed to version 2.0.
//
// Revision 1.4  2007/09/01 23:40:14  uid6503
// Corrected:
// 1) compilation issues
// 	lib/covfit/proplist.cc
// 	main/utils/merge_proplists.cc
// 2) error in source-sink parsing
// 	main/strippers/strip_hadspec.cc
//
// Revision 1.3  2007/04/16 21:31:31  kostas
// fixed fitter automake
//
// Revision 1.2  2007/04/10 14:16:45  kostas
// mege_proplists now merges into complex numbers
//
// Revision 1.1  2007/03/22 05:11:39  kostas
// added the merge_proplist utility
//
// merges proplist files into a single multicolumn file
//

#include <iostream>
#include <cstdio>

#include "covfit/proplist.h"

using namespace std;
using namespace CovFit ;

int main(int argc, char **argv)
{
  if ((argc < 3))
  {
    cerr << "Usage: " << argv[0] << " <out proplist>  <proplist> <proplist> ...\n" ;
    cerr << "   or  " << argv[0] << " -c <out proplist>  <re proplist> <im proplist> ...\n" ;
    cerr << "   It merges proplist files into a single multicolumn file\n" ;
    cerr << "   If the -c flag is specified then it takes pairs of real and imaginary part proplists and then merges them into complex numbers in a multicolumn file\n" ;
    cerr << endl;
    exit(1);
  }

  int Nfiles = (argc - 2);
  int first_file = 2 ;
  int incr  = 1;
  if(string(argv[1]) == "-c"){
    cout<<"Merging into comlex numbers\n" ;
    Nfiles = Nfiles - 1;
    first_file = 3 ;
    incr = 2 ;
    if(Nfiles%2){
      cerr<<"Need even number of files for complex merging\n";
      exit(2);
    }
  }


  cout << "Number of files: "<<Nfiles<<endl ;
  PropList pp ;
  Array< Array< Array<Double> > > av ;
  Array< Array< Array<DComplex> > > c ;

  int Nx,Nt ;
  for(int i(first_file);i<argc;i++){
    int ttNx=ReadProplist(pp,argv[i]) ;
    int ttNt = pp[0].size() ;
    if(i==first_file){
      Nx=ttNx ;
      Nt=ttNt ;
      //cout<<"av.size()="<<av.size()<<endl;
      cout<<"pp.size()="<<pp.size()<<endl;
      if(incr==1){
	av.resize(pp.size());
	for(int j(0);j<pp.size();j++){
	  av[j].resize(pp[j].size());
	  for(int k(0);k<av[j].size();k++){
	    av[j][k].resize(Nfiles) ;
	  }
	}
	cout<<"av.size()="<<av.size()<<endl;
	cout<<"av[0].size()="<<av[0].size()<<endl;
	cout<<"av[0][0].size()="<<av[0][0].size()<<endl;
      }
      else{
	c.resize(pp.size());
	for(int j(0);j<pp.size();j++){
	  c[j].resize(pp[j].size());
	  for(int k(0);k<c[j].size();k++){
	    c[j][k].resize(Nfiles/2) ;
	  }
	}
	cout<<"c.size()="<<c.size()<<endl;
	cout<<"c[0].size()="<<c[0].size()<<endl;
	cout<<"c[0][0].size()="<<c[0][0].size()<<endl;
      }
    }

    if(Nx!=ttNx){
      cerr<<"Nx: "<<Nx<<" ttNx: "<<ttNx<<endl ;
      cerr<<"OOPS! wrong proplist header\n" ;
      exit(2);
    }
    if(Nt!=ttNt){ 
      cerr<<"Nt: "<<Nt<<" ttNt: "<<ttNt<<endl ;
      cerr<<"OOPS! wrong propagator length\n" ;
      exit(1);
    }
    
    if(incr == 1){
      for(int j(0);j<pp.size();j++)
	for(int k(0);k<av[j].size();k++)
	  av[j][k][i-2]=pp[j][k];
    }
    else{
      //cout<<i<<" OK\n";
      if( i%2 == 1)
	for(int j(0);j<pp.size();j++)
	  for(int k(0);k<c[j].size();k++)
	    c[j][k][(i-3)/2]=DComplex( pp[j][k], imag(c[j][k][(i-3)/2]) );
      if( i%2 == 0)
	for(int j(0);j<pp.size();j++)
	  for(int k(0);k<c[j].size();k++)
	    c[j][k][(i-3)/2]=DComplex( real(c[j][k][(i-3)/2]), pp[j][k] );
    }
  
  }


  if(incr == 1) 
    WriteProplist(av,Nx,argv[first_file-1]);
  else
    WriteProplist(c ,Nx,argv[first_file-1]);

  /** DEBUG **
  for(int j(0);j<c.size();j++)
    for(int k(0);k<c[0].size();k++){
      cout<<k<<" ";
      for(int i(0);i<c[0][0].size();i++)
	cout<<c[j][k][i].re<<" "<<c[j][k][i].im<<" ";
      cout<<endl ;
    }
  **/
}

