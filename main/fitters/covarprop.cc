//$Id: covarprop.cc,v 2.0 2008/12/05 04:43:46 edwards Exp $
//$Log: covarprop.cc,v $
//Revision 2.0  2008/12/05 04:43:46  edwards
//Changed to version 2.0.
//
//Revision 1.2  2005/04/27 19:53:52  edwards
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
//Revision 1.1  2004/02/13 15:32:18  kostas
//Added Doug's fitting programs
//

/*! simple dirver pogram that computes and writes to disk
   the covariance matrix
 */


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "covfit/covarmat.h"

using namespace CovFit;
using namespace std;

int main(int argc,char *argv[]){
  int blocksize,npoints ;
 
  /* find size of propagator */
  if(argc < 3){fprintf(stderr,"Not enough arguments\n"); exit(1);}
  if(sscanf(argv[1],"%d",&npoints) != 1){
    fprintf(stderr,"Bad argument (maximum distance)\n"); exit(1);
  }
  if(sscanf(argv[2],"%d",&blocksize) != 1){
    fprintf(stderr,"Bad argument (blocksize)\n"); exit(1);
  }
  npoints++;	/* zero distance counts! */

  CovarMat cm(npoints,blocksize);
  cm.ReadDataList("stdin");
  cm.CalcCovarMat();
  cm.WriteCovarMat("stdout");
  return 0;
}
