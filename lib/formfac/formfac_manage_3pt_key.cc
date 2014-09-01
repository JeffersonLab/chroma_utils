// $Id: formfac_manage_3pt_key.cc,v 2.1 2009/03/21 21:33:45 edwards Exp $
/*! \file
 * \brief Key for manage 3-pt funcs
 */

#include "formfac/formfac_manage_3pt_key.h"

#include "adat/map_traits.h"

#undef FF_DEBUG

namespace FF
{
  using namespace ADAT;

  //-----------------------------------------------------------------------------------
  //! Reader of a threept arg
  void read(XMLReader& xml, const std::string& path, ThreePtArg& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "PiPf", param.pi_pf);
    read(paramtop, "gamma", param.gamma);
    read(paramtop, "links", param.links);
    read(paramtop, "src_name", param.src_name);
    read(paramtop, "src_smear", param.src_smear);
    read(paramtop, "src_lorentz", param.src_lorentz);
    read(paramtop, "src_spin", param.src_spin);
    read(paramtop, "snk_name", param.snk_name);
    read(paramtop, "snk_smear", param.snk_smear);
    read(paramtop, "snk_lorentz", param.snk_lorentz);
    read(paramtop, "snk_spin", param.snk_spin);
    read(paramtop, "dt", param.dt);
    read(paramtop, "quark", param.quark);
    read(paramtop, "mass", param.mass);
    read(paramtop, "ensemble", param.ensemble);
  }

  //! Writer of a threept arg
  void write(XMLWriter& xml, const std::string& path, const ThreePtArg& param)
  {
    push(xml, path);

    write(xml, "PiPf", param.pi_pf);
    write(xml, "gamma", param.gamma);
    write(xml, "links", param.links);
    write(xml, "src_name", param.src_name);
    write(xml, "src_smear", param.src_smear);
    write(xml, "src_lorentz", param.src_lorentz);
    write(xml, "src_spin", param.src_spin);
    write(xml, "snk_name", param.snk_name);
    write(xml, "snk_smear", param.snk_smear);
    write(xml, "snk_lorentz", param.snk_lorentz);
    write(xml, "snk_spin", param.snk_spin);
    write(xml, "quark", param.quark);
    write(xml, "dt", param.dt);
    write(xml, "mass", param.mass);
    write(xml, "ensemble", param.ensemble);

    pop(xml);
  }


  //-----------------------------------------------------------------------------------
  //! Reader of a threept arg
  void read(BinaryReader& bin, ThreePtArg& param)
  {
    read(bin, param.pi_pf.p_i);
    read(bin, param.pi_pf.p_f);
    read(bin, param.gamma);
    read(bin, param.links);
    readDesc(bin, param.src_name);
    readDesc(bin, param.src_smear);
    read(bin, param.src_lorentz);
    read(bin, param.src_spin);
    readDesc(bin, param.snk_name);
    readDesc(bin, param.snk_smear);
    read(bin, param.snk_lorentz);
    read(bin, param.snk_spin);
    read(bin, param.dt);
    read(bin, param.quark);
    readDesc(bin, param.mass);
    readDesc(bin, param.ensemble);
  }

  //! Writer of a threept arg
  void write(BinaryWriter& bin, const ThreePtArg& param)
  {
    write(bin, param.pi_pf.p_i);
    write(bin, param.pi_pf.p_f);
    write(bin, param.gamma);
    write(bin, param.links);
    writeDesc(bin, param.src_name);
    writeDesc(bin, param.src_smear);
    write(bin, param.src_lorentz);
    write(bin, param.src_spin);
    writeDesc(bin, param.snk_name);
    writeDesc(bin, param.snk_smear);
    write(bin, param.snk_lorentz);
    write(bin, param.snk_spin);
    write(bin, param.quark);
    write(bin, param.dt);
    writeDesc(bin, param.mass);
    writeDesc(bin, param.ensemble);
  }

  //-----------------------------------------------------------------------------------
  // Merge this key with another key to produce a final key
  void mergeKeys(ThreePtArg& out, const ThreePtArgReduced& in, const ThreePtArg& prototype)
  {
    out.pi_pf       = in.pi_pf;
    out.gamma       = in.gamma;
    out.links       = in.links;

    out.dt          = prototype.dt;
    out.quark       = prototype.quark;

    out.src_name    = prototype.src_name;
    out.src_smear   = prototype.src_smear;
    out.src_lorentz = prototype.src_lorentz;
    out.src_spin    = prototype.src_spin;

    out.snk_name    = prototype.snk_name;
    out.snk_smear   = prototype.snk_smear;
    out.snk_lorentz = prototype.snk_lorentz;
    out.snk_spin    = prototype.snk_spin;

    out.mass        = prototype.mass;
    out.ensemble    = prototype.ensemble;
  }


} // namespace FF
