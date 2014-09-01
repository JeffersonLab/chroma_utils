// $Id: formfac_solver_row.cc,v 2.0 2008/12/05 04:43:37 edwards Exp $
/*! \file
 * \brief Keys used for solver
 */

#include "formfac/formfac_solver_row.h"

namespace FF 
{
  // Namespace composition
  using namespace ENSEM;

  // Read parameters
  void read(XMLReader& xml, const std::string& path, TimeFitRange_t& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "t_i", param.t_i);
    read(paramtop, "t_f", param.t_f);
  }

  // Writer
  void write(XMLWriter& xml, const std::string& path, const TimeFitRange_t& param)
  {
    push(xml, path);

    write(xml, "t_i", param.t_i);
    write(xml, "t_f", param.t_f);

    pop(xml);
  }

} // namespace FF
