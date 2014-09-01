//$Id: adat_util.h,v 2.0 2008/12/05 04:43:37 edwards Exp $
//$Log: adat_util.h,v $
//Revision 2.0  2008/12/05 04:43:37  edwards
//Changed to version 2.0.
//
//Revision 1.1  2008/11/12 21:37:21  edwards
//Big change for the library. Moved away from "include" and "lib" model.
//Now, everthing is in a "lib" directory. Really not a big change for
//main programs. The corresponding main/*/Makefile.am are pointed to "lib"
//instead of "include".
//
//Revision 1.1  2005/07/29 14:35:16  edwards
//Changed include of files now down in "io" subdir.
//
//Revision 1.4  2005/04/27 19:53:51  edwards
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
//Revision 1.3  2005/02/11 03:11:49  kostas
//added a block data utility
//
//Revision 1.2  2004/02/10 22:21:05  kostas
//Added some statistics routines an mres computation and some operators
//for the Array class
//
// -*- C++ -*-
//
#ifndef __ADAT__UTILS_H___
#define __ADAT__UTILS_H___

namespace ADATUtil
{
  //! Is the native byte order big endian?
  bool big_endian();

  //! Byte-swap an array of data each of size nmemb
  void byte_swap(void *ptr, size_t size, size_t nmemb);

  //! fread on a binary file written in big-endian order
  size_t bfread(void *ptr, size_t size, size_t nmemb, FILE *stream);

  //! fwrite to a binary file in big-endian order
  size_t bfwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);
}

#endif
