// $Id: strip_hadron_contract.cc,v 2.0 2008/12/05 04:44:03 edwards Exp $

#include <list>
#include "hadron_contract_factory.h"
#include "hadron_contract_aggregate.h"
#include "io/adat_qdpio.h"

//using namespace ADAT;
//using namespace Util;
using namespace Strippers;

// Linkage hack
bool linkage()
{
  bool success = true;
  // Register all factories
  success &= Strippers::HadronContractEnv::registerAll();

  return success;
}


//
// Main program - loop over files
//
int main(int argc, char **argv)
{
  // Register all factories
  linkage();
  
  // Grab arguments
  if (argc <= 1)
  {
    std::cerr << "Usage: " << argv[0] << " file1 [file2 file3... fileN]" << std::endl;
    exit(1);
  }

  std::string record_root = "/HadronContraction";
  std::string record_id   = "/ContractionType";

  // The number of data files (size of the ensemble)
  int nbin = argc - 1;
  
  std::cout << "Number measurements = " << nbin << std::endl;

  // The list of measurements - this is a list for each record
  std::list < Handle< HadronContract > > the_measurements;

  // Read the first file to setup all the object factories
  try 
  {
    std::cout << "Open file " << argv[1] << std::endl;
    XMLReader file_xml;
    ADATFileReader qio_in(file_xml, argv[1]);

    int nrec = 0;
    while(qio_in.eof())
    {
      // Read the record
      std::cout << "Read record " << nrec << std::endl;
      XMLReader xml;
      Handle<HadronContractResult_t> had_cont;
      qio_in.read(xml, had_cont->bin);
      {
	std::ostringstream os;
	xml.print(os);
	had_cont->xml = os.str();
      }

      // Extract the id for this record
      std::string id;
      read(xml, record_id, id);

      // Construct the factory for this id
      Handle<HadronContract> hadronContract(
	TheHadronContractFactory::Instance().createObject(
	  id,
	  xml,
	  record_root,
	  nbin));

      // Since we have the binary, process it
      (*hadronContract)(had_cont, 0);

      // Put this object on the list. This will copy the handle.
      the_measurements.push_back(hadronContract);
      ++nrec;
    }
  }
  catch(std::bad_cast) 
  {
    std::cerr << "CHROMA: caught cast error" << std::endl;
    exit(1);
  }
  catch(std::bad_alloc) 
  { 
    // This might happen on any node, so report it
    std::cerr << argv[0] << ": caught bad memory allocation" << std::endl;
    exit(1);
  }
  catch(const std::string& e) 
  {
    std::cerr << argv[0] << ": Caught Exception: " << e << std::endl;
    exit(1);
  }
  catch(std::exception& e) 
  {
    std::cerr << argv[0] << ": Caught standard library exception: " << e.what() << std::endl;
    exit(1);
  }
  catch(...)
  {
    // This might happen on any node, so report it
    std::cerr << argv[0] << ": caught generic exception during measurement" << std::endl;
    exit(1);
  }


  // Read and process the rest of the files
  for(int ibin=2; ibin < argc; ++ibin)
  {
    const std::string filename = argv[ibin+1];
    try 
    {
      std::cout << "Open file " << filename << std::endl;
      XMLReader file_xml;
      ADATFileReader qio_in(file_xml, filename);

      // Read the records
      int nrec = 0;
      for(std::list< Handle<HadronContract> >::const_iterator had_ptr= the_measurements.begin(); 
	  had_ptr != the_measurements.end(); 
	  ++had_ptr, ++nrec)
      {
	if (qio_in.eof())
	{
	  std::cerr << argv[0] << ": unexpected end-of-file in file= " << filename << std::endl;
	  exit(1);
	}

	Handle<HadronContract> hadronContract = *had_ptr;

	// Read the record
	std::cout << "Read record " << nrec << std::endl;
	Handle<HadronContractResult_t> had_cont;

	{
	  XMLReader xml;
	  qio_in.read(xml, had_cont->bin);

	  std::ostringstream os;
	  xml.print(os);
	  had_cont->xml = os.str();
	}
	
	// Since we have the binary, process it
	(*hadronContract)(had_cont, ibin);
      }
    }
    catch(std::bad_cast) 
    {
      std::cerr << "CHROMA: caught cast error" << std::endl;
      exit(1);
    }
    catch(std::bad_alloc) 
    { 
      // This might happen on any node, so report it
      std::cerr << argv[0] << ": caught bad memory allocation" << std::endl;
      exit(1);
    }
    catch(const std::string& e) 
    {
      std::cerr << argv[0] << ": Caught Exception: " << e << std::endl;
      exit(1);
    }
    catch(std::exception& e) 
    {
      std::cerr << argv[0] << ": Caught standard library exception: " << e.what() << std::endl;
      exit(1);
    }
    catch(...)
    {
      // This might happen on any node, so report it
      std::cerr << argv[0] << ": caught generic exception during measurement" << std::endl;
      exit(1);
    }
  } // for ibin

  // Write the files
// print_files(files, spec, false);

  exit(0);
}
