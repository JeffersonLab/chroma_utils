//
/*! @file
 * @brief XML IO support
 */

#include "ensem/ensem.h"

using std::iostream;
using std::string;

namespace ENSEM {

//--------------------------------------------------------------------------------
//! Read a XML Array element
template<typename T>
void readArrayPrimitive(XMLReader& xml, const std::string& s, Array<T>& result)
{
  std::ostringstream error_message;
  
  // Try reading the list as a string
  string list_string;
  read(xml, s, list_string);

  // Count the number of elements
  std::istringstream list_stream(list_string);
	
  int array_size = 0;
  T dummy;
  while(list_stream >> dummy)
    ++array_size;

  if ((! list_stream.eof()) && list_stream.fail())
  {
    error_message << "Error in reading array " << s << std::endl;
    throw error_message.str();
  }

  if (array_size == 0)
  {
    error_message << "Something wrong with reading array " << list_string << std::endl;
    throw error_message.str();
  }
      
  // Now resize the array to hold the no of elements.
  result.resize(array_size);

  // Get the elements one by one
  // I do not understand why, but use a new stringstream
//  list_stream.str(list_string);
  std::istringstream list_stream2(list_string);

  for(int i=0; i < result.size(); i++) 
  {
    // read the element.
    list_stream2 >> result[i];
  }
}

void read(XMLReader& xml, const std::string& xpath, Array<Integer>& result)
{
  readArrayPrimitive<Integer>(xml, xpath, result);
}
void read(XMLReader& xml, const std::string& xpath, Array<Real32>& result)
{
  readArrayPrimitive<Real32>(xml, xpath, result);
}
void read(XMLReader& xml, const std::string& xpath, Array<Real64>& result)
{
  readArrayPrimitive<Real64>(xml, xpath, result);
}
void read(XMLReader& xml, const std::string& xpath, Array<Boolean>& result)
{
  readArrayPrimitive<Boolean>(xml, xpath, result);
}


//--------------------------------------------------------------------------------
// Write an array of basic types
template<typename T>
void writeArrayPrimitive(XMLWriter& xml, const std::string& s, const Array<T>& s1)
{
  std::ostringstream output;

  output << s1[0];
  for(int index=1; index < s1.size(); index++) 
    output << " " << s1[index];
    
  // Write the array - do not use a normal string write
  xml.openTag(s);
  xml << output.str();
  xml.closeTag();
}


void write(XMLWriter& xml, const std::string& xpath, const Array<Integer>& output)
{
  writeArrayPrimitive<Integer>(xml, xpath, output);
}
void write(XMLWriter& xml, const std::string& xpath, const Array<Real32>& output)
{
  writeArrayPrimitive<Real32>(xml, xpath, output);
}
void write(XMLWriter& xml, const std::string& xpath, const Array<Real64>& output)
{
  writeArrayPrimitive<Real64>(xml, xpath, output);
}
void write(XMLWriter& xml, const std::string& xpath, const Array<Boolean>& output)
{
  writeArrayPrimitive<Boolean>(xml, xpath, output);
}

}  // namespace ENSEM
