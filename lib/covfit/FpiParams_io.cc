//$Id: FpiParams_io.cc,v 2.0 2008/12/05 04:43:32 edwards Exp $
//$Log: FpiParams_io.cc,v $
//Revision 2.0  2008/12/05 04:43:32  edwards
//Changed to version 2.0.
//
//Revision 1.2  2005/05/13 16:46:49  kostas
//added fpi
//
//Revision 1.1  2005/05/13 16:38:17  kostas
//added fpi stuff
//
//


#include "covfit/FpiParams_io.h"

using namespace std;
 
namespace CovFit {
  void read(XMLReader& xml, const string& path, FpiParams_t& fp)
  {
    XMLReader top(xml, path);
    
    read(top,"Fitter"  ,fp.fit);
    read(top,"CovarMat",fp.cov);
    read(top,"filenameSP",fp.filenameSP);
    read(top,"filenameSS",fp.filenameSS);
    
    read(top,"quark_mass",fp.quark_mass);
    
    read(top,"mint",fp.mint);
    read(top,"maxt",fp.maxt);
    
    read(top,"fit_paramsSP",fp.fit_paramsSP);
    read(top,"fit_paramsSS",fp.fit_paramsSS);
    
    read(top,"mres",fp.mres);
    read(top,"block",fp.block);
  }
  
  
  void write(XMLWriter& xml, const string& path, const FpiParams_t& fp)
  {
    
    write(xml,"Fitter"  ,fp.fit);
    write(xml,"CovarMat",fp.cov);
    write(xml,"filenameSP",fp.filenameSP);
    write(xml,"filenameSS",fp.filenameSS);
    
    write(xml,"quark_mass",fp.quark_mass);
    
    write(xml,"mint",fp.mint);
    write(xml,"maxt",fp.maxt);
    
    write(xml,"fit_paramsSP",fp.fit_paramsSP);
    write(xml,"fit_paramsSS",fp.fit_paramsSS);
    
    write(xml,"mres",fp.mres);
    write(xml,"block",fp.block);
  }

}
