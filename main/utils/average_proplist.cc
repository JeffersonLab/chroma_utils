// $Id: average_proplist.cc,v 2.0 2008/12/05 04:44:06 edwards Exp $
// $Log: average_proplist.cc,v $
// Revision 2.0  2008/12/05 04:44:06  edwards
// Changed to version 2.0.
//
// Revision 1.2  2005/10/30 20:10:03  kostas
// fixed average_proplist
//
// Revision 1.1  2005/10/30 19:58:20  kostas
// adde average_proplist
//
//

#include <iostream>
#include <cstdio>

#include "covfit/proplist.h"

using namespace std;
using namespace CovFit ;

int main(int argc, char **argv)
{
  if ((argc < 4)||(argc%2==1))
  {
    cerr << "Usage: " << argv[0] << " <out proplist> <sign> <proplist> <sign> <proplist> ...\n" ;
    cerr << "   It computes (<sign> <proplist> <sign> <proplist> ...)/<number of files>\n" ;
    cerr << endl;
    exit(1);
  }

  int Nfiles = (argc - 1)/2;
  cout << "Number of files: "<<Nfiles<<endl ;
  PropList pp ;
  PropList av ;

  int Nx,Nt ;
  for(int i(2);i<argc;i+=2){
    int ttNx=ReadProplist(pp,argv[i+1]) ;
    int ttNt = pp[0].size() ;
    if(i==2){
      Nx=ttNx ;
      Nt=ttNt ;
      cout<<"av.size()="<<av.size()<<endl;
      cout<<"pp.size()="<<pp.size()<<endl;
      av.resize(pp.size());
      for(int j(0);j<pp.size();j++){
	av[j].resize(pp[j].size());
	av[j] = 0.0 ;
      }
      cout<<"av.size()="<<av.size()<<endl;
      cout<<"av[0].size()="<<av[0].size()<<endl;
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

    
    if(argv[i] == string("+") )
      av += pp ;
    else if( argv[i] == string("-") )
      av -= pp ;
    else{
      cerr<<"OOPS! unknown sign |"<< argv[i]<<"|"<<endl ;
      exit(1) ;
    }
  }

  for(int j(0);j<av.size();j++){
    av[j] /= Real(Nfiles);
  }
	     

  WriteProplist(av,Nx,argv[1]);
  

}

