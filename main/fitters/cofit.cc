//$Id: cofit.cc,v 2.0 2008/12/05 04:43:46 edwards Exp $
//$Log: cofit.cc,v $
//Revision 2.0  2008/12/05 04:43:46  edwards
//Changed to version 2.0.
//
//Revision 1.5  2005/11/23 04:23:03  kostas
//...
//
//Revision 1.4  2005/11/23 03:28:55  kostas
//fixed the path to fitters in cofit.cc
//
//Revision 1.3  2005/04/27 19:53:52  edwards
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
//Revision 1.2  2004/03/17 23:34:44  kostas
//small fix bugs
//
//Revision 1.1  2004/02/13 15:32:18  kostas
//Added Doug's fitting programs
//

/*! simple dirver pogram that does covariant fits */


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
  
  CovarMat cv;
  cv.ReadCovarMat(argv[1]);
  coFitter fit(P,para,cv);
  fit.Fit();
  fit.graph("cofit_graph.ax");

}
