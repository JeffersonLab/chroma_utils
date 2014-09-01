// $Id: formfac_manage_E_key.cc,v 2.1 2009/03/21 21:33:45 edwards Exp $
/*! \file
 * \brief Keys for manage energy funcs
 */

#include "adat/map_traits.h"
#include "formfac/formfac_manage_E_key.h"

namespace FF
{
  using namespace ADAT;

  //-----------------------------------------------------------------------------------
  //! Reader of a energy arg
  void read(XMLReader& xml, const std::string& path, EnergyArg& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "level", param.level);
    read(paramtop, "mom", param.mom);
    read(paramtop, "mass", param.mass);
    read(paramtop, "ensemble", param.ensemble);
  }

  //! Writer of a energy arg
  void write(XMLWriter& xml, const std::string& path, const EnergyArg& param)
  {
    push(xml, path);

    write(xml, "level", param.level);
    write(xml, "mom", param.mom);
    write(xml, "mass", param.mass);
    write(xml, "ensemble", param.ensemble);

    pop(xml);
  }


  //-----------------------------------------------------------------------------------
  //! Reader of a energy arg
  void read(BinaryReader& bin, EnergyArg& param)
  {
    read(bin, param.level);
    read(bin, param.mom);
    readDesc(bin, param.mass);
    readDesc(bin, param.ensemble);
  }

  //! Writer of a energy arg
  void write(BinaryWriter& bin, const EnergyArg& param)
  {
    write(bin, param.level);
    write(bin, param.mom);
    writeDesc(bin, param.mass);
    writeDesc(bin, param.ensemble);
  }

} // namespace FF
