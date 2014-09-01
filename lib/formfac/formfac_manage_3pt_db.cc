// $Id: formfac_manage_3pt_db.cc,v 2.2 2009/03/05 04:23:50 edwards Exp $
//
// Manage building-blocks


#include "formfac/formfac_manage_3pt_db.h"
#include "formfac/formfac_manage_factory.h"

namespace FF
{
  //----------------------------------------------------------------------------------
  // Parameters
  Manage3PtFuncDBParams_t::Manage3PtFuncDBParams_t(XMLReader& xml, const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "db_file", db_file);
    read(paramtop, "PrototypeKey", proto_key);
    read(paramtop, "LatticeParam", lattice);
  }


  // Parameters
  void Manage3PtFuncDBParams_t::writeXML(XMLWriter& xml, const std::string& path) const
  {
    push(xml, path);

    write(xml, "db_file", db_file);
    write(xml, "PrototypeKey", proto_key);
    write(xml, "LatticeParam", lattice);

    pop(xml);
  }


  // Parameters
  void read(XMLReader& xml, const std::string& path, Manage3PtFuncDBParams_t& param)
  {
    Manage3PtFuncDBParams_t tmp(xml, path);
    param = tmp;
  }


  // Parameters
  void write(XMLWriter& xml, const std::string& path, const Manage3PtFuncDBParams_t& param)
  {
    param.writeXML(xml, path);
  }


  //----------------------------------------------------------------------------------
  //! Constructor
  Manage3PtFuncDB::Manage3PtFuncDB(const std::string& dbase, 
				   const KeyHadron3PtCorr_t& proto_key_,
				   const LatticeParam& lattice_) :
    proto_key(proto_key_),
    lattice(lattice_)
  {
    if (database.open(dbase, O_RDONLY, 0644) != 0)
    {
      std::cerr << __func__ << ": error opening dbase= " << dbase << std::endl;
      exit(1);
    }
  }


  //! Destructor
  Manage3PtFuncDB::~Manage3PtFuncDB() {}

  //! Return a value given a key
  EnsemVectorComplexF Manage3PtFuncDB::operator[](const ThreePtArg& arg)
  {
    // Convert the input argument key with the prototype key
    K ky;
    mergeKeys(ky, arg, proto_key);

    SerialDBKey<K> key;
    key.key() = ky;

    std::vector< SerialDBData<SV> > vals;
    int ret;
    if ((ret = database.get(key, vals)) != 0)
    {
      std::cerr << __func__ << ": key not found\n" << ky;
      exit(1);
    }

    V eval;
    eval.resize(vals.size());
    eval.resizeObs(vals[0].data().numElem());

    for(int i=0; i < vals.size(); ++i)
    {
      SV sval = vals[i].data();
      pokeEnsem(eval, sval, i);
    }

    return eval;
  }

  //! Number of bins
  int Manage3PtFuncDB::size() const
  {
    return database.size();
  }

  //! Time extent
  int Manage3PtFuncDB::timeLen() const
  {
    return lattice.latt_size[lattice.decay_dir];
  }

  //! Decay direction
  int Manage3PtFuncDB::decayDir() const
  {
    return lattice.decay_dir;
  }

  //! Lattice size
  const ArrayInt& Manage3PtFuncDB::lattSize() const
  {
    return lattice.latt_size;
  }


  //! Manager factory
  namespace FormfacManage3PtFuncDBEnv
  { 
    //! Anonymous namespace
    namespace
    {
      bool registered = false;

      const std::string name("THREEPT_DB");

      //-------------------- callback functions ---------------------------------------

      //! Callback
      Manage3PtFunc* createFunc(XMLReader& xml_in,
				const std::string& path)
      {
	Manage3PtFuncDBParams_t params(xml_in, path);

	return new Manage3PtFuncDB(params.db_file, 
				   params.proto_key,
				   params.lattice);
      }
    }  // anonymous namespace


    // Register the callbacks
    bool registerAll(void) 
    {
      bool success = true;

      if (! registered)
      {
	//! Register all the factories
	success &= TheManage3PtFuncFactory::Instance().registerObject(name, createFunc);

	registered = true;
      }

      return success;
    }

  }  // end namespace FormfacManageEnv

} // namespace FF
