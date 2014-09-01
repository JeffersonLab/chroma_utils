// $Id: make_mom_list.cc,v 2.0 2008/12/05 04:43:48 edwards Exp $

#include "ensem/ensem.h"
#include "formfac/formfac_ensemble.h"
#include "formfac/formfac_qsq.h"

using namespace ENSEM;
using namespace FF;
using namespace std;

int main(int argc, char *argv[])
{
  double r_1_latt = 2.662      ; // in lattice units
  double r_1_phys = 0.324      ; // in fm

  LatticeParam lattice;
  lattice.latt_size.resize(Nd);
  lattice.latt_size         = 12;
  lattice.latt_size[Nd-1]   = 48;
  lattice.decay_dir    = Nd-1;
  lattice.a_fm         =  0.10;  // fm
  lattice.aniso.anisoP = true;
  lattice.aniso.c_sq   = 1.04;
  lattice.aniso.xi     = 3;
  
  QsqParam param;
  param.nq_norm_sq_max = 25;
  param.lattice        = lattice;
  param.mass_i         = 0.485;
  param.mass_f         = 0.485;

  // Available sink momentum
  {
    ArrayInt mm(3);
    mm[0] =  0; mm[1] =  0; mm[2] =  0; param.nf_list.push_back(mm);
    mm[0] =  1; mm[1] =  0; mm[2] =  0; param.nf_list.push_back(mm);
    mm[0] = -1; mm[1] = -1; mm[2] =  0; param.nf_list.push_back(mm);
//    mm[0] = -1; mm[1] = -1; mm[2] = -1; param.nf_list.push_back(mm);
//    mm[0] = -2; mm[1] =  0; mm[2] =  0; param.nf_list.push_back(mm);
  }

  // Construct map
  std::list<QsqVal> qsq_val = constructQsqFixedMom(param);

  // Print out map
  cout << qsq_val;
}

  
