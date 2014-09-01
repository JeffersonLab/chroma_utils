// $Id: formfac_ensemble.cc,v 2.0 2008/12/05 04:43:35 edwards Exp $

#include <map>
#include "formfac/formfac_qsq.h"

//using namespace ENSEM;
//using namespace std;

namespace FF
{

  // Read in params
  void read(XMLReader& xml, const std::string& path, AnisoParam& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "anisoP", param.anisoP);
    read(paramtop, "c_sq", param.c_sq);
    read(paramtop, "xi", param.xi);
    read(paramtop, "xi_0", param.xi_0);
    read(paramtop, "nu", param.nu);
    read(paramtop, "m0_at", param.m0_at);
  }

  // Write params
  void write(XMLWriter& xml, const std::string& path, const AnisoParam& param)
  {
    push(xml, path);

    write(xml, "anisoP", param.anisoP);
    write(xml, "c_sq", param.c_sq);
    write(xml, "xi", param.xi);
    write(xml, "xi_0", param.xi_0);
    write(xml, "nu", param.nu);
    write(xml, "m0_at", param.m0_at);

    pop(xml);
  }


  // Read in params
  void read(XMLReader& xml, const std::string& path, LatticeParam& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "latt_size", param.latt_size);
    read(paramtop, "decay_dir", param.decay_dir);
    read(paramtop, "a_fm", param.a_fm);

    if (paramtop.count("AnisoParam") != 0)
    {
      read(paramtop, "AnisoParam", param.aniso);
    }

    // Lots of sanity checks
    if (param.latt_size.size() == 0)
    {
      std::cerr << __func__ << ": invalid lattice_size" << std::endl;
      exit(1);
    }

    for(int i=0; i < param.latt_size.size(); ++i)
    {
      if (param.latt_size[i] <= 0)
      {
	std::cerr << __func__ << ": invalid lattice_size extent = " << param.latt_size[i] << std::endl;
	exit(1);
      }
    }

    if (param.decay_dir < 0 || param.decay_dir > param.latt_size.size())
    {
      std::cerr << __func__ << ": invalid decay direction" << std::endl;
      exit(1);
    }
  }

  // Writer params
  void write(XMLWriter& xml, const std::string& path, const LatticeParam& param)
  {
    push(xml, path);

    write(xml, "latt_size", param.latt_size);
    write(xml, "decay_dir", param.decay_dir);
    write(xml, "a_fm", param.a_fm);

    if (param.aniso.anisoP)
      write(xml, "AnisoParam", param.aniso);

    pop(xml);
  }


  //! Reader
  void read(XMLReader& xml, const std::string& path, PiPf& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "p_f", param.p_f);
    read(paramtop, "p_i", param.p_i);
  }

  //! Writer
  void write(XMLWriter& xml, const std::string& path, const PiPf& param)
  {
    push(xml, path);

    write(xml, "p_f", param.p_f);
    write(xml, "p_i", param.p_i);
    write(xml, "q_sq", norm2(param.p_f - param.p_i));

    pop(xml);
  }


  //! Reader
  void read(XMLReader& xml, const std::string& path, EnsembleInfo& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "nbin", param.nbin);
    read(paramtop, "LatticeParam", param.lattice);

  }

  //! Writer
  void write(XMLWriter& xml, const std::string& path, const EnsembleInfo& param)
  {
    push(xml, path);

    write(xml, "nbin", param.nbin);
    write(xml, "LatticeParam", param.lattice);

    pop(xml);
  }


  //! Support for maps
  bool operator<(const PiPf& a, const PiPf& b)
  {
    return concat(a.p_i, a.p_f) < concat(b.p_i, b.p_f);
  }

} // namespace FF
