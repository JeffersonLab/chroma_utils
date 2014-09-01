// -*- C++ -*-
// $Id: current_ops.h,v 2.0 2008/12/05 04:43:47 edwards Exp $
/*! \file
 * \brief Current operators
 */

#ifndef __current_ops_h__
#define __current_ops_h__

#include "adat/handle.h"
#include "io/adat_xml_group_reader.h"
#include "formfac/formfac_ensemble.h"
#include "formfac/formfac_manage_3pt.h"
#include "formfac/formfac_manage_factory.h"
#include "formfac/formfac_solver_row.h"

namespace FF
{
  using namespace ADAT;
  using namespace ENSEM;

  //----------------------------------------------------------------------------------
  //! Base class for Minkowski-space current operators
  class CurrentOperator
  {
  public:  
    //! Virtual destructor
    virtual ~CurrentOperator() {}

    //! The current operator combination of 3-pt correlators
    virtual EnsemReal operator()(const LLSqRow_t& row) const = 0;

  protected:
    //! Initialize a 3-pt manage object
    virtual Manage3PtFunc* initialize3Pt(const GroupXML_t& manage_xml) const;
  };


  //----------------------------------------------------------------------------------
  //! Params for insertion operator
  struct CurrentOperatorParams
  {
    CurrentOperatorParams();
    CurrentOperatorParams(XMLReader& in, const std::string& path);
    void writeXML(XMLWriter& in, const std::string& path) const;

    GroupXML_t         threept;            /*!< xml holding 3-pt manage object */
    LatticeParam       lattice;            /*!< Holds lattice size and aniso*/
  };


  //! Reader
  void read(XMLReader& xml, const std::string& path, CurrentOperatorParams& param);

  //! Writer
  void write(XMLWriter& xml, const std::string& path, const CurrentOperatorParams& param);


  //----------------------------------------------------------------------------------
  //! Local vector current
  /*!
   * Matrix element is  V_\mu
   */
  class LocalVectorCurrentOperator : public CurrentOperator
  {
  public:
    //! Full constructor
    LocalVectorCurrentOperator(const CurrentOperatorParams& p);

    //! Destructor
    ~LocalVectorCurrentOperator() {}

    //! The current operator combination of 3-pt correlators
    EnsemReal operator()(const LLSqRow_t& row) const;

  private:
    CurrentOperatorParams  params;   /*!< params */
    Handle<Manage3PtFunc>  threept;
  };



  //----------------------------------------------------------------------------------
  //! Conserved vector current
  /*!
   * Matrix element is  V_\mu
   */
  class ConservedVectorCurrentOperator : public CurrentOperator
  {
  public:
    //! Full constructor
    ConservedVectorCurrentOperator(const CurrentOperatorParams& p);

    //! Destructor
    ~ConservedVectorCurrentOperator() {}

    //! The current operator combination of 3-pt correlators
    EnsemReal operator()(const LLSqRow_t& row) const;

  private:
    CurrentOperatorParams  params;   /*!< params */
    Handle<Manage3PtFunc>  threept;
  };



  //----------------------------------------------------------------------------------
  //! q_upol_0
  class QUpol0CurrentOperator : public CurrentOperator
  {
  public:
    //! Full constructor
    QUpol0CurrentOperator(const CurrentOperatorParams& p);

    //! Destructor
    ~QUpol0CurrentOperator() {}

    //! The current operator combination of 3-pt correlators
    EnsemReal operator()(const LLSqRow_t& row) const;

  private:
    CurrentOperatorParams  params;   /*!< params */
    Handle<Manage3PtFunc>  threept;
  };



} // namespace FF

#endif
