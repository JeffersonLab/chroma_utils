// $Id: ensem_io.cc,v 2.0 2008/12/05 04:43:33 edwards Exp $
//
// QDP data parallel interface
//

#include <stdlib.h>
#include "ensem/ensem_io.h"

using std::iostream;
using std::string;

namespace ENSEM {

//-----------------------------------------
//! text reader support
TextReader::TextReader() {}

TextReader::TextReader(const string& p) {open(p);}

void TextReader::open(const string& p) 
{
  f.open(p.c_str(),std::ifstream::in);

  if (! is_open())
  {
    std::cerr << "TextReader: failed to open file: " << p << std::endl;
    exit(1);
  }
}

void TextReader::close()
{
  if (is_open()) 
    f.close();
}

// Propagate status to all nodes
bool TextReader::is_open()
{
  return f.is_open();
}

// Propagate status to all nodes
bool TextReader::fail()
{
  return f.fail();
}

TextReader::~TextReader() {close();}


// String reader
void TextReader::read(string& input)
{
  getIstream() >> input;
}

// Readers
void TextReader::read(char& input) 
{
  readPrimitive<char>(input);
}
void TextReader::read(int& input) 
{
  readPrimitive<int>(input);
}
void TextReader::read(unsigned int& input)
{
  readPrimitive<unsigned int>(input);
}
void TextReader::read(short int& input)
{
  readPrimitive<short int>(input);
}
void TextReader::read(unsigned short int& input)
{
  readPrimitive<unsigned short int>(input);
}
void TextReader::read(long int& input)
{
  readPrimitive<long int>(input);
}
void TextReader::read(unsigned long int& input)
{
  readPrimitive<unsigned long int>(input);
}
void TextReader::read(float& input)
{
  readPrimitive<float>(input);
}
void TextReader::read(double& input)
{
  readPrimitive<double>(input);
}
void TextReader::read(bool& input)
{
  readPrimitive<bool>(input);
}

template< typename T>
void TextReader::readPrimitive(T& input)
{
  getIstream() >> input;
}

// Different bindings for read functions
TextReader& operator>>(TextReader& txt, string& input)
{
  txt.read(input);
  return txt;
}
TextReader& operator>>(TextReader& txt, char& input)
{
  txt.read(input);
  return txt;
}
TextReader& operator>>(TextReader& txt, int& input)
{
  txt.read(input);
  return txt;
}
TextReader& operator>>(TextReader& txt, unsigned int& input)
{
  txt.read(input);
  return txt;
}
TextReader& operator>>(TextReader& txt, short int& input)
{
  txt.read(input);
  return txt;
}
TextReader& operator>>(TextReader& txt, unsigned short int& input)
{
  txt.read(input);
  return txt;
}
TextReader& operator>>(TextReader& txt, long int& input)
{
  txt.read(input);
  return txt;
}
TextReader& operator>>(TextReader& txt, unsigned long int& input)
{
  txt.read(input);
  return txt;
}
TextReader& operator>>(TextReader& txt, float& input)
{
  txt.read(input);
  return txt;
}
TextReader& operator>>(TextReader& txt, double& input)
{
  txt.read(input);
  return txt;
}
TextReader& operator>>(TextReader& txt, bool& input)
{
  txt.read(input);
  return txt;
}


//-----------------------------------------
//! text writer support
TextWriter::TextWriter() {}

TextWriter::TextWriter(const string& p) {open(p);}

void TextWriter::open(const string& p)
{
  f.open(p.c_str(),std::ofstream::out | std::ofstream::trunc);

  if (! is_open())
  {
    std::cerr << "TextWriter: failed to open file: " << p << std::endl;
    exit(1);
  }
}

void TextWriter::close()
{
  if (is_open()) 
    f.close();
}

// Propagate status to all nodes
bool TextWriter::is_open()
{
  return f.is_open();
}

void TextWriter::flush()
{
  f.flush();
}

// Propagate status to all nodes
bool TextWriter::fail()
{
  return f.fail();
}

TextWriter::~TextWriter() {close();}


// Different bindings for write functions
TextWriter& operator<<(TextWriter& txt, const string& output)
{
  txt.write(output);
  return txt;
}

TextWriter& operator<<(TextWriter& txt, const char* output)
{
  txt.write(output);
  return txt;
}

TextWriter& operator<<(TextWriter& txt, int output)
{
  txt.write(output);
  return txt;
}

TextWriter& operator<<(TextWriter& txt, unsigned int output)
{
  txt.write(output);
  return txt;
}

TextWriter& operator<<(TextWriter& txt, short int output)
{
  txt.write(output);
  return txt;
}

TextWriter& operator<<(TextWriter& txt, unsigned short int output)
{
  txt.write(output);
  return txt;
}

TextWriter& operator<<(TextWriter& txt, long int output)
{
  txt.write(output);
  return txt;
}

TextWriter& operator<<(TextWriter& txt, unsigned long int output)
{
  txt.write(output);
  return txt;
}

TextWriter& operator<<(TextWriter& txt, float output)
{
  txt.write(output);
  return txt;
}

TextWriter& operator<<(TextWriter& txt, double output)
{
  txt.write(output);
  return txt;
}

TextWriter& operator<<(TextWriter& txt, bool output)
{
  txt.write(output);
  return txt;
}




void TextWriter::write(const string& output)
{
  writePrimitive<string>(output);
}

void TextWriter::write(const char* output)
{
  write(string(output));
}

void TextWriter::write(const char& output) 
{
  writePrimitive<char>(output);
}

void TextWriter::write(const int& output) 
{
  writePrimitive<int>(output);
}

void TextWriter::write(const unsigned int& output)
{
  writePrimitive<unsigned int>(output);
}

void TextWriter::write(const short int& output)
{
  writePrimitive<short int>(output);
}

void TextWriter::write(const unsigned short int& output)
{
  writePrimitive<unsigned short int>(output);
}

void TextWriter::write(const long int& output)
{
  writePrimitive<long int>(output);
}

void TextWriter::write(const unsigned long int& output)
{
  writePrimitive<unsigned long int>(output);
}

void TextWriter::write(const float& output)
{
  writePrimitive<float>(output);
}

void TextWriter::write(const double& output)
{
  writePrimitive<double>(output);
}

void TextWriter::write(const bool& output)
{
  writePrimitive<bool>(output);
}

template< typename T>
void TextWriter::writePrimitive(const T& output)
{
  getOstream() << output;
}



} // namespace ENSEM
