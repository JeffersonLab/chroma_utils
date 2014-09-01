//$Id: fit.cc,v 2.0 2008/12/05 04:43:46 edwards Exp $
//$Log: fit.cc,v $
//Revision 2.0  2008/12/05 04:43:46  edwards
//Changed to version 2.0.
//
//Revision 1.6  2005/11/23 04:23:03  kostas
//...
//
//Revision 1.5  2005/11/23 03:18:51  kostas
//moved the fitter header files into the main include tree
//
//Revision 1.4  2005/04/27 19:53:52  edwards
//Completely overhauled include and lib structure. Put all
//of Kostas' fit stuff into their own covfit subdir in include
//and lib. All in its own namespace. Similarly, all of
//ensem is in it's own dirs and namespaces. Unified use of
//XML Array as base array type. I've not touched TTT which
//uses someother classes.
//The new convention should be is people put there little mini
//packages in their own dirs and namespaces. One can compose
//a bigger namespace with namespace composition.
//
//Revision 1.3  2004/03/17 22:30:51  kostas
//modified how to compile fitfunctions
//
//Revision 1.2  2004/02/26 20:37:56  kostas
//fixed problems with automake
//
//Revision 1.1  2004/02/13 15:32:18  kostas
//Added Doug's fitting programs
//

/*! simple dirver pogram that does fits */

#include <cstdio>
#include <iostream>
#include <string>
#include <covfit/fitters/fitfunctions.h>
#include <covfit/fitter.h>

using namespace CovFit;
using namespace std;

int main(int argc,char *argv[])
{
  FUNCTION P;
  if( argc!=2  &&  argc!=3) {
    printf("\tfit error\n\tcorrect format is:\n");
    printf("\tfit   datafile   [parameterfile]\n");
    exit(0);
  }
  char data[200];
  char para[200];
  strcpy(para,"STDIN") ;
  if(argc>=2){
    strcpy(data,argv[1]) ;
  }
  if(argc==3){
    strcpy(para,argv[2]) ;
  }
  
  Fitter fit(P,para,data);
  fit.Fit();
  fit.graph("fit_graph.ax");
}
