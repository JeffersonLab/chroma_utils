//$Id: readcovp.cc,v 2.0 2008/12/05 04:43:46 edwards Exp $
//$Log: readcovp.cc,v $
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

/*! simple dirver pogram that reads the covatiance matrix from
  a file and writes it to the stdout
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "covfit/covarmat.h"

using namespace CovFit;
using namespace std;

int main(int argc,char *argv[])
{
  CovarMat cm;
  cm.ReadCovarMat(argv[1]);
  cm.CalcCovarMat();
  cm.WriteCovarMat("stdout");
  return 0;
}
