//  $Id: proginfo.cc,v 1.2 2009/04/23 03:28:47 edwards Exp $
/*! \file
 *  \brief Print out basic info about this program
 */

#include "io/adat_xmlio.h"
#include <ctime>


namespace ADAT
{
  using namespace ADATXML;

  //! Print out basic information about this program
  /*!
   *  \param xml          The XML stream to which to write the information.
   */

  void proginfo(XMLWriter& xml)
  {
    push(xml,"ProgramInfo");

    push(xml,"code_version");
    write(xml, "chroma", PACKAGE_STRING);
    pop(xml);

    std::time_t now;

    if(std::time(&now)==-1)
    {
      std::cerr << __func__ << ": Cannot get the time.\n";
      return;
    }
    std::tm *tp = std::localtime(&now);

    char date[64];
    std::strftime(date, 63, "%d %b %y %X %Z", tp);
    write(xml,"run_date", date);

    pop(xml);
  }

}  // end namespace Chroma
