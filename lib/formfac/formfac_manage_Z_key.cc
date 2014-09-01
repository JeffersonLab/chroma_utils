// $Id: formfac_manage_Z_key.cc,v 2.1 2009/03/21 21:33:45 edwards Exp $
/*! \file
 * \brief Manage amplitudes
 */

#include "adat/map_traits.h"
#include "formfac/formfac_manage_Z_key.h"

namespace FF
{
  using namespace ADAT;
 
  //-----------------------------------------------------------------------------------
  //! Reader of a amp arg
  void read(XMLReader& xml, const std::string& path, AmpArg& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "level", param.level);
    read(paramtop, "mom", param.mom);
    read(paramtop, "name", param.name);
    read(paramtop, "smear", param.smear);
    read(paramtop, "mass", param.mass);
    read(paramtop, "ensemble", param.ensemble);
  }

  //! Writer of a amp arg
  void write(XMLWriter& xml, const std::string& path, const AmpArg& param)
  {
    push(xml, path);

    write(xml, "level", param.level);
    write(xml, "mom", param.mom);
    write(xml, "name", param.name);
    write(xml, "smear", param.smear);
    write(xml, "mass", param.mass);
    write(xml, "ensemble", param.ensemble);

    pop(xml);
  }


  //-----------------------------------------------------------------------------------
  //! Reader of a amp arg
  void read(BinaryReader& bin, AmpArg& param)
  {
    read(bin, param.level);
    read(bin, param.mom);
    readDesc(bin, param.name);
    readDesc(bin, param.smear);
    readDesc(bin, param.mass);
    readDesc(bin, param.ensemble);
  }

  //! Writer of a amp arg
  void write(BinaryWriter& bin, const AmpArg& param)
  {
    write(bin, param.level);
    write(bin, param.mom);
    writeDesc(bin, param.name);
    writeDesc(bin, param.smear);
    writeDesc(bin, param.mass);
    writeDesc(bin, param.ensemble);
  }

} // namespace FF
