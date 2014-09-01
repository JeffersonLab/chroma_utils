// $Id: merge_cproplists.cc,v 2.0 2008/12/05 04:44:06 edwards Exp $
// $Log: merge_cproplists.cc,v $
// Revision 2.0  2008/12/05 04:44:06  edwards
// Changed to version 2.0.
//
// Revision 1.1  2007/05/10 20:42:42  kostas
// added merge
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
  if ((argc < 2))
  {
    cerr << "Usage: " << argv[0] << " <out proplist>  <proplist> <proplist> ...\n" ;
    cerr << "   It merges proplist files into a single multicolumn file\n" ;
    cerr << "   If the -c flag is specified then it takes pairs of real and imaginary part proplists and then merges them into complex numbers in a multicolumn file\n" ;
    cerr << endl;
    exit(1);
  }

  int Nfiles = (argc - 2);
  int first_file = 2 ;

  cout << "Number of files: "<<Nfiles<<endl ;
  cPropList pp ;
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
      c.resize(pp.size());
      for(int j(0);j<pp.size();j++){
	c[j].resize(pp[j].size());
	for(int k(0);k<c[j].size();k++){
	  c[j][k].resize(Nfiles) ;
	}
      }
      cout<<"c.size()="<<c.size()<<endl;
      cout<<"c[0].size()="<<c[0].size()<<endl;
      cout<<"c[0][0].size()="<<c[0][0].size()<<endl;
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
    
    for(int j(0);j<pp.size();j++)
      for(int k(0);k<c[j].size();k++)
	c[j][k][i-2]=pp[j][k];
 
  }

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

