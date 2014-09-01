/*! \file
 *  \brief Read an XML group as a string
 */

#include "io/adat_xml_group_reader.h"
#include <vector>

namespace ADATXML
{

  // Read group and return in a struct
  GroupXML_t readXMLGroup(XMLReader& xml_in, 
			  const std::string& path, const std::string& type_name)
  {
    GroupXML_t group;

    try
    {
      XMLReader xml_tmp(xml_in, path);
      std::ostringstream os;
      xml_tmp.print(os);
      read(xml_tmp, type_name, group.id);
      group.xml = os.str();
      group.path = "/" + path;
    }
    catch(const std::string& e) 
    {
      std::cerr << __func__ << ": caught exception reading XML: " << e << std::endl;
      exit(1);
    }

    return group;
  }


  // Read group and return in a struct
  Array<GroupXML_t> readXMLArrayGroup(XMLReader& xml_in, 
				      const std::string& path, 
				      const std::string& type_name)
  {
    Array<GroupXML_t> group;

    try
    {
      XMLReader xml_tmp(xml_in, path);
      group.resize(xml_tmp.count("elem"));

      for(int i=0; i < group.size(); i++) 
      {
	// Create the query for the element 
	std::ostringstream element_xpath;
	element_xpath << "elem[" << (i+1) << "]";

	XMLReader xml_elem(xml_tmp, element_xpath.str());
	std::ostringstream os;
	xml_elem.print(os);
	read(xml_elem, type_name, group[i].id);
	group[i].xml = os.str();
	group[i].path = "/elem";
      }
    }
    catch(const std::string& e) 
    {
      std::cerr << __func__ << ": caught exception reading XML: " << e << std::endl;
      exit(1);
    }

    return group;
  }

  // Read group and return in a struct
  std::vector<GroupXML_t> readXMLVectorGroup(XMLReader& xml_in, 
					     const std::string& path, 
					     const std::string& type_name)
  {
    std::vector<GroupXML_t> group;

    try
    {
      XMLReader xml_tmp(xml_in, path);
      group.resize(xml_tmp.count("elem"));

      for(int i=0; i < group.size(); i++) 
      {
	// Create the query for the element 
	std::ostringstream element_xpath;
	element_xpath << "elem[" << (i+1) << "]";

	XMLReader xml_elem(xml_tmp, element_xpath.str());
	std::ostringstream os;
	xml_elem.print(os);
	read(xml_elem, type_name, group[i].id);
	group[i].xml = os.str();
	group[i].path = "/elem";
      }
    }
    catch(const std::string& e) 
    {
      std::cerr << __func__ << ": caught exception reading XML: " << e << std::endl;
      exit(1);
    }

    return group;
  }

}  // end namespace Chroma
