// $Id: thin_bb.cc,v 2.0 2008/12/05 04:43:49 edwards Exp $
// Thin out a BB file

#include "formfac/formfac_bb.h"
#include <stdlib.h>

using namespace FF;
using namespace std;

int main(int argc, char *argv[])
{
  if (argc != 4)
  {
    cerr << "Usage:"
	 << argv[0] << ":  <desired links_max>  <input BB>  <output BB>" 
	 << endl;
    exit(1);
  }

  int new_links_max = atoi(argv[1]);

  // Read the old BB file
  BuildingBlocks_t old_bar;
  read(argv[2], old_bar);
  
  // Thin the puppy
  BuildingBlocks_t new_bar;
  thinBB(new_bar, new_links_max, old_bar);

  // Now write
  write(argv[3], new_bar);
  
//  cout << argv[0] << ": exiting" << endl;
  exit(0);
}

