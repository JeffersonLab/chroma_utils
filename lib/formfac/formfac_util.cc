// $Id: formfac_util.cc,v 2.0 2008/12/05 04:43:37 edwards Exp $
//
// Utilities

#include "ensem/ensem.h"
#include "formfac/formfac_util.h"

namespace FF
{

  //! Compute energy
  EnsemReal energyDisp(const EnsemReal& m, const ArrayInt& p, const LatticeParam& lattice)
  {
    const double c_sq      =  lattice.aniso.c_sq;
    const double xi        =  lattice.aniso.xi;
    const double xi_sq     =  xi * xi;

    Array<double> nnf = contMom(p,lattice.latt_size[0]);
    Array<Real> nf(nnf.size());
    for(int i=0; i < nf.size(); ++i)
      nf[i] = nnf[i];

    EnsemReal  e = sqrt( Real(c_sq/xi_sq)*norm2(nf) + m*m);

    return e;
  }

  //! special photon energy routine
  EnsemReal photonenergyDisp(const EnsemReal& Qsq, const ArrayInt& p, const LatticeParam& lattice)
  {
    //const double c_sq      =  lattice.aniso.c_sq;
    // !don't use the hadron dipersion relation??????
    const double xi        =  lattice.aniso.xi;
    const double xi_sq     =  xi * xi;

    Array<double> nnf = contMom(p,lattice.latt_size[0]);
    Array<Real> nf(nnf.size());
    for(int i=0; i < nf.size(); ++i)
      nf[i] = nnf[i];


    EnsemReal  e = sqrt( Real(1.0/xi_sq)*norm2(nf) - Qsq);
 

    return e;
  }

  //! Decompose a lexicographic site into coordinates
  Array<int> crtesn(int ipos, const Array<int>& latt_size)
  {
    Array<int> coord(latt_size.size());

    /* Calculate the Cartesian coordinates of the VALUE of IPOS where the 
     * value is defined by
     *
     *     for i = 0 to NDIM-1  {
     *        X_i  <- mod( IPOS, L(i) )
     *        IPOS <- int( IPOS / L(i) )
     *     }
     *
     * NOTE: here the coord(i) and IPOS have their origin at 0. 
     */
    for(int i=0; i < latt_size.size(); ++i)
    {
      coord[i] = ipos % latt_size[i];
      ipos = ipos / latt_size[i];
    }

    return coord;
  }


  //! Calculates the lexicographic site index from the coordinate of a site
  /*! 
   * Nothing specific about the actual lattice size, can be used for 
   * any kind of latt size 
   */
  int local_site(const Array<int>& coord, const Array<int>& latt_size)
  {
    int order = 0;

    for(int mmu=latt_size.size()-1; mmu >= 1; --mmu)
      order = latt_size[mmu-1]*(coord[mmu] + order);

    order += coord[0];

    return order;
  }


} // namespace FF
