//$Id: FitParams_io.cc,v 2.0 2008/12/05 04:43:32 edwards Exp $
//$Log: FitParams_io.cc,v $
//Revision 2.0  2008/12/05 04:43:32  edwards
//Changed to version 2.0.
//
//Revision 1.4  2007/03/03 04:02:43  kostas
//fixed spectrum to accept periodic functions
//
//Revision 1.3  2006/12/08 06:48:35  kostas
//fixed bug in write state
//
//Revision 1.2  2006/12/08 05:23:06  kostas
//modified spectrum to admit arbitrary function
//
//Revision 1.1  2005/04/27 19:53:52  edwards
//Completely overhauled include and lib structure. Put all
//of Kostas' fit stuff into their own covfit subdir in include
//and lib. All in its own namespace. Similarly, all of
//ensem is in it's own dirs and namespaces. Unified use of
//XML Array as base array type. I've not touched TTT which
//uses someother classes.
//The new convention should be is people put there little mini
//packages in their own dirs and namespaces. One can compose
//a bigger namespace with namespace composition.
//
//Revision 1.3  2004/10/18 18:16:22  kostas
//converted the particle database to an xmlarray
//
//Revision 1.2  2004/05/20 16:33:35  kostas
//cvs tags
//

#include "covfit/FitParams_io.h"

using namespace std;
 
namespace CovFit {

  void read(XMLReader& xml, const std::string& path, FitterParam& c){
    XMLReader top(xml, path);
    
    read(top,"fitfunc",c.fitfunc);
    
    read(top,"min_dist",c.min_dist);
    read(top,"max_dist",c.max_dist);
    
    read(top,"fit_params",c.fit_params);
    read(top,"mass_param",c.mass_param);
    if(top.count("fixedpar")!=0)
      read(top,"fixedpar",c.fixedpar);

    c.period = -1969.02 ;
    if(top.count("period")!=0)
      read(top,"period",c.period);
  }
  
  void write(XMLWriter& top, const std::string& path, const FitterParam& c){
    push(top,path);
    
    write(top,"fitfunc",c.fitfunc);
    write(top,"min_dist",c.min_dist);
    write(top,"max_dist",c.max_dist);
    
    write(top,"fit_params",c.fit_params);
    write(top,"mass_param",c.mass_param);
    
    
    if(c.fixedpar.size()!=0)
      write(top,"fixedpar",c.fixedpar);
    
    if(c.period!= -1969.02)
      write(top,"period",c.period);

    pop(top);
    
  }


  void read(XMLReader& xml, const string& path, FitParams_t& fp)
  {
    XMLReader top(xml, path);
    
    read(top,"Fitter"  ,fp.fit);
    read(top,"CovarMat",fp.cov);
    read(top,"States"  ,fp.sta);   
  }

  void read(XMLReader& xml, const string& path, DBFitParams_t& fp)
  {
    XMLReader top(xml, path);
    
    read(top,"Fitter"  ,fp.fit);
    read(top,"CovarMat",fp.cov);
    read(top,"States"  ,fp.sta);   
  }
  
  
  void read(XMLReader& xml, const string& path, Fitter_t& ft)
  {
    XMLReader top(xml, path);
    
    read(top,"Niter"  ,ft.Niter);
    read(top,"Toler"  ,ft.Toler);
  }
  
  void read(XMLReader& xml, const string& path, CovarMat_t& cv)
  {
    XMLReader top(xml, path);
    
    read(top,"block"  ,cv.block);
    read(top,"fold"   ,cv.fold);
  }
  
  
  //! Read a sink type enum
  void read(XMLReader& xml, const string& path, SinkType& param)
  {
    string src_type_str;
    read(xml, path, src_type_str);
    if (src_type_str == "POINT_SINK")
      param = SNK_TYPE_POINT_SINK;
    else if (src_type_str == "WALL_SINK")
      param = SNK_TYPE_WALL_SINK;
    else if (src_type_str == "SHELL_SINK")
      param = SNK_TYPE_SHELL_SINK;
    else if (src_type_str == "BNDST_SINK")
      param = SNK_TYPE_BNDST_SINK;
    else 
      {
	cerr << "Unsupported SinkType" << endl;
	exit(1);
      }
  }


  //! Read a source type enum
  void read(XMLReader& xml, const string& path, SourceType& param)
  {
    string src_type_str;
    read(xml, path, src_type_str);
    if (src_type_str == "POINT_SOURCE")
      param = SRC_TYPE_POINT_SOURCE;
    else if (src_type_str == "WALL_SOURCE")
      param = SRC_TYPE_WALL_SOURCE;
    else if (src_type_str == "SHELL_SOURCE")
      param = SRC_TYPE_SHELL_SOURCE;
    else if (src_type_str == "BNDST_SOURCE")
      param = SRC_TYPE_BNDST_SOURCE;
    else 
      {
	cerr << "Unsupported SourceType" << endl;
	exit(1);
      }
  }

  
  void read(XMLReader& xml, const string& path, State_t& st)
  {
    XMLReader top(xml, path);
    
    read(top,"name"   ,st.name);
    read(top,"filename"   ,st.filename);
    
    read(top,"source" ,st.source);
    read(top,"sink"   ,st.sink);
    
    read(top,"quark_mass",st.quark_mass);
    
    if(st.quark_mass.size()==2) st.type = MESON ;
    if(st.quark_mass.size()==3) st.type = BARYON ;
    
    read(top,"stopper",st.stopper);
    read(top,"range",st.range);
    
    read(top,"fit",st.fit);
    
  }
  

  //! Write a source type enum
  void write(XMLWriter& xml, const string& path, SourceType param)
  {
    string src_type_str;
    if (param == SRC_TYPE_POINT_SOURCE)
      src_type_str = "POINT_SOURCE";
    else if (param == SRC_TYPE_WALL_SOURCE)
      src_type_str = "WALL_SOURCE";
    else if (param == SRC_TYPE_SHELL_SOURCE)
      src_type_str = "SHELL_SOURCE";
    else if (param == SRC_TYPE_BNDST_SOURCE)
      src_type_str = "BNDST_SOURCE";
    else 
      {
	cerr << "Unsupported SourceType" << endl;
	exit(1);
      }
    write(xml, path, src_type_str);
  }

  void read(XMLReader& xml, const string& path, DBState_t& st)
  {
    XMLReader top(xml, path);
    
    read(top,"name"   ,st.name);
    read(top,"filename"   ,st.filename);
        
    read(top,"stopper",st.stopper);
    read(top,"range",st.range);

    read(top,"Key",st.key);
    
    read(top,"fit",st.fit);
    
  }
  

  
  
  //! Write a sink type enum
  void write(XMLWriter& xml, const string& path, SinkType param)
  {
    string src_type_str;
    if (param == SNK_TYPE_POINT_SINK)
      src_type_str = "POINT_SINK";
    else if (param == SNK_TYPE_WALL_SINK)
      src_type_str = "WALL_SINK";
    else if (param == SNK_TYPE_SHELL_SINK)
      src_type_str = "SHELL_SINK";
    else if (param == SNK_TYPE_BNDST_SINK)
      src_type_str = "BNDST_SINK";
    else 
      {
	cerr << "Unsupported SinkType" << endl;
	exit(1);
      }
    write(xml, path, src_type_str);
  }
  
  
  
  
  void write(XMLWriter& xml, const string& path, const FitParams_t& fp)
  {
    push(xml, path);
    
    write(xml,"Fitter"  ,fp.fit);
    write(xml,"CovarMat",fp.cov);
    write(xml,"States"   ,fp.sta);   
    
    pop(xml) ;
  }

  void write(XMLWriter& xml, const string& path, const DBFitParams_t& fp)
  {
    push(xml, path);
    
    write(xml,"Fitter"  ,fp.fit);
    write(xml,"CovarMat",fp.cov);
    write(xml,"States"   ,fp.sta);   
    
    pop(xml) ;
  }

  void write(XMLWriter& xml, const string& path, const Fitter_t& ft)
  {
    push(xml, path);
    
    write(xml, "Niter", ft.Niter);
    write(xml, "Toler", ft.Toler);
    
    pop(xml) ;
  }
  
  void write(XMLWriter& xml, const string& path, const CovarMat_t& cv)
  {
    push(xml, path);
    
    write(xml,"block"  ,cv.block);
    write(xml,"fold"   ,cv.fold);
    
    pop(xml) ;
  }
  
  void write(XMLWriter& xml, const string& path, const State_t& st)
  {
    
    push(xml, path);
    
    write(xml,"name"     ,st.name); 
    write(xml,"filename" ,st.filename);
    write(xml,"source"   ,st.source);
    write(xml,"sink"     ,st.sink);
    
    write(xml,"quark_mass",st.quark_mass);
    write(xml,"stopper",st.stopper);
    write(xml,"range",st.range);
    
    write(xml,"fit",st.fit);
   
    
    pop(xml) ;
    
  }

  void write(XMLWriter& xml, const string& path, const DBState_t& st)
  {
    
    push(xml, path);
    
    write(xml,"name"     ,st.name); 
    write(xml,"filename" ,st.filename);
    write(xml,"key" ,st.key);
    
    write(xml,"stopper",st.stopper);
    write(xml,"range",st.range);
    
    write(xml,"fit",st.fit);
   
    
    pop(xml) ;
    
  }

} // namespace CovFit
