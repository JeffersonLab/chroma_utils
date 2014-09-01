// $Id: formfac_manage_2pt_key.cc,v 2.2 2009/03/26 21:21:43 edwards Exp $
//
// Manage building-blocks


#include "formfac/formfac_manage_2pt.h"
#include "formfac/formfac_qsq.h"

namespace FF
{
  //-----------------------------------------------------------------------------------
  // Reader of a twopt arg
  void read(XMLReader& xml, const std::string& path, TwoPtArg& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "mom", param.mom);
    read(paramtop, "src_name", param.src_name);
    read(paramtop, "src_smear", param.src_smear);
    read(paramtop, "src_lorentz", param.src_lorentz);
    read(paramtop, "src_spin", param.src_spin);
    read(paramtop, "snk_name", param.snk_name);
    read(paramtop, "snk_smear", param.snk_smear);
    read(paramtop, "snk_lorentz", param.snk_lorentz);
    read(paramtop, "snk_spin", param.snk_spin);
    read(paramtop, "mass", param.mass);
    read(paramtop, "ensemble", param.ensemble);
  }

  // Writer of a twopt arg
  void write(XMLWriter& xml, const std::string& path, const TwoPtArg& param)
  {
    push(xml, path);

    write(xml, "mom", param.mom);
    write(xml, "src_name", param.src_name);
    write(xml, "src_smear", param.src_smear);
    write(xml, "src_lorentz", param.src_lorentz);
    write(xml, "src_spin", param.src_spin);
    write(xml, "snk_name", param.snk_name);
    write(xml, "snk_smear", param.snk_smear);
    write(xml, "snk_lorentz", param.snk_lorentz);
    write(xml, "snk_spin", param.snk_spin);
    write(xml, "mass", param.mass);
    write(xml, "ensemble", param.ensemble);

    pop(xml);
  }


  //-----------------------------------------------------------------------------------
  // Reader of a twopt arg
  void read(BinaryReader& bin, TwoPtArg& param)
  {
    read(bin, param.mom);
    readDesc(bin, param.src_name);
    readDesc(bin, param.src_smear);
    read(bin, param.src_lorentz);
    read(bin, param.src_spin);
    readDesc(bin, param.snk_name);
    readDesc(bin, param.snk_smear);
    read(bin, param.snk_lorentz);
    read(bin, param.snk_spin);
    readDesc(bin, param.mass);
    readDesc(bin, param.ensemble);
  }

  // Writer of a twopt arg
  void write(BinaryWriter& bin, const TwoPtArg& param)
  {
    write(bin, param.mom);
    writeDesc(bin, param.src_name);
    writeDesc(bin, param.src_smear);
    write(bin, param.src_lorentz);
    write(bin, param.src_spin);
    writeDesc(bin, param.snk_name);
    writeDesc(bin, param.snk_smear);
    write(bin, param.snk_lorentz);
    write(bin, param.snk_spin);
    writeDesc(bin, param.mass);
    writeDesc(bin, param.ensemble);
  }

} // namespace FF
