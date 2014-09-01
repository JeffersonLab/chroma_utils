#
# Path and flags for all files used in building the library
#
AM_CXXFLAGS = -I@top_srcdir@/lib \
              -I@top_srcdir@/other_libs/c-lime/include \
              -I@top_srcdir@/other_libs/xpath_reader/include \
              -I@top_builddir@/other_libs/c-lime/include \
              -I@top_srcdir@/other_libs/filedb/src \
              -I@top_srcdir@/other_libs/filedb/filehash \
              -I@top_srcdir@/include/recipes

AM_CXXFLAGS += @CXXFLAGS@ @LIBXML2_CXXFLAGS@

#
# Include files
#
#nodist_include_HEADERS =

# Adat
nobase_include_HEADERS = adat/handle.h \
    adat/funcmap.h  adat/objfactory.h  adat/objfunctor.h  \
    adat/typeinfo.h  adat/typemanip.h  adat/singleton.h  \
    adat/typelist.h  adat/typetraits.h \
    adat/adat_stopwatch.h  adat/map_traits.h  adat/map_obj.h \
    adat/proginfo.h adat/adat_arrays.h adat/adat_arraynd.h

# Io
nobase_include_HEADERS += io/adat_xmlio.h io/adat_io.h \
    io/adat_byteorder.h io/adat_xml_group_reader.h io/key_val_db.h \
    io/old_key_val_db.h

# Covfit
nobase_include_HEADERS += covfit/statistics.h covfit/proplist.h \
    covfit/Function.h covfit/fitter.h covfit/covarmat.h \
    covfit/jackknife.h \
    covfit/FitParams_io.h covfit/FpiParams_io.h \
    covfit/fitters/chiPTfp.h        covfit/fitters/exponper.h \
    covfit/fitters/exponaper_bar.h  covfit/fitters/fitfunctions.h \
    covfit/fitters/exponaper.h      covfit/fitters/polynomium.h \
    covfit/fitters/expon.h     covfit/fitters/Lookup.h \
    covfit/fitters/expon_altern.h covfit/fitters/exponper_altern.h \
    covfit/fitters/exponper_altern_onemass.h \
    covfit/fitters/chiPTlog.h \
    covfit/fitters/BaryonFfunc.h covfit/fitters/BaryonSfunc.h \
    covfit/fitters/gAJfunc.h covfit/fitters/gAKfunc.h

# Ensemble
nobase_include_HEADERS += ensem/ensem.h ensem/ensem_ensem.h \
    ensem/ensem_forward.h ensem/ensem_io.h \
    ensem/ensem_newops.h ensem/ensem_observable.h \
    ensem/ensem_obsmatrix.h ensem/ensem_obsscalar.h \
    ensem/ensem_obstensor.h ensem/ensem_obsvector.h \
    ensem/ensem_operator_tags.h ensem/ensem_params.h \
    ensem/ensem_precision.h ensem/ensem_primgamma.h \
    ensem/ensem_primitive.h ensem/ensem_primmatrix.h \
    ensem/ensem_primscalar.h ensem/ensem_primseed.h \
    ensem/ensem_primspinmat.h ensem/ensem_primspinvec.h \
    ensem/ensem_primcolormat.h ensem/ensem_primcolorvec.h \
    ensem/ensem_primvector.h ensem/ensem_random.h \
    ensem/ensem_reality.h ensem/ensem_scalar_specific.h \
    ensem/ensem_scalarsite_defs.h ensem/ensem_simpleword.h \
    ensem/ensem_specializations.h \
    ensem/ensem_traits.h ensem/ensem_type_computations.h

# Numerical recipes
nobase_include_HEADERS += \
    recipes/nr.h recipes/nrtypes.h recipes/nrtypes_nr.h \
    recipes/nrutil.h recipes/nrutil_nr.h recipes/print_array.h

# Formfactor support
nobase_include_HEADERS += \
    formfac/formfac.h formfac/formfac_ensemble.h \
    formfac/formfac_solver_row.h \
    formfac/formfac_solver.h formfac/formfac_solver_driver.h \
    formfac/formfac_qsq.h \
    formfac/formfac_util.h formfac/formfac_bb.h formfac/formfac_state.h \
    formfac/formfac_manage.h \
    formfac/formfac_manage_2pt.h \
    formfac/formfac_manage_2pt_key.h \
    formfac/formfac_manage_2pt_stripped.h \
    formfac/formfac_manage_E.h \
    formfac/formfac_manage_E_key.h \
    formfac/formfac_manage_E_stripped.h \
    formfac/formfac_manage_Z.h \
    formfac/formfac_manage_Z_key.h \
    formfac/formfac_manage_Z_stripped.h \
    formfac/formfac_manage_npr.h \
    formfac/formfac_manage_3pt.h \
    formfac/formfac_manage_3pt_key.h \
    formfac/formfac_manage_3pt_db.h \
    formfac/formfac_manage_aggregate.h \
    formfac/formfac_manage_factory.h \
    formfac/formfac_meson_fit.h \
    formfac/formfac_data_covar.h \
    formfac/formfac_ensem_linfit.h \
    formfac/formfac_linear_covar_genfit.h \
    formfac/formfac_chisqq.h \
    formfac/hadron_1pt_corr.h \
    formfac/hadron_2pt_corr.h \
    formfac/hadron_3pt_corr.h

# Parton distribution sources
nobase_include_HEADERS += \
    parton/parton_correlation_functions.h parton/parton_derivatives.h \
    parton/parton_distribution_moments.h parton/parton_renormalization_factors.h


#
# Library
#
lib_LIBRARIES = libadat.a

# Adat simple utility sources
libadat_a_SOURCES = adat/adat_stopwatch.cc adat/map_traits.cc \
    adat/proginfo.cc

# Io sources
libadat_a_SOURCES += io/adat_xmlio.cc io/adat_io.cc io/adat_byteorder.cc \
    io/adat_crc32.cc io/adat_xml_group_reader.cc

# CovFit sources
libadat_a_SOURCES += covfit/proplist.cc covfit/coFitter.cc \
    covfit/CovarMat.cc covfit/Fitter.cc covfit/BaseFitter.cc \
    covfit/FitParams_io.cc  covfit/FpiParams_io.cc \
    covfit/function_switchboard.cc

# Ensemble source
libadat_a_SOURCES += ensem/ensem_random.cc ensem/ensem_io.cc ensem/ensem_xmlio.cc

# Numerical recipes
libadat_a_SOURCES += \
    recipes/brent.cc recipes/spline.cc recipes/splint.cc \
    recipes/gammln.cc recipes/gammp.cc recipes/gammq.cc \
    recipes/gcf.cc recipes/gser.cc \
    recipes/pythag.cc \
    recipes/svbksb.cc recipes/svdcmp.cc recipes/svdfit.cc recipes/svdvar.cc \
    recipes/dawson.cc recipes/gauleg.cc recipes/factrl.cc recipes/plgndr.cc \
    recipes/zbrak.cc recipes/zbrac.cc recipes/zbrent.cc recipes/mnbrak.cc \
    recipes/mnewt.cc recipes/ludcmp.cc recipes/lubksb.cc

# Formfac sources
libadat_a_SOURCES += \
    formfac/formfac_ensemble.cc \
    formfac/formfac_solver_row.cc \
    formfac/formfac_solver.cc formfac/formfac_solver_driver.cc \
    formfac/formfac_qsq.cc \
    formfac/formfac_util.cc formfac/formfac_bb.cc \
    formfac/formfac_manage_2pt_key.cc \
    formfac/formfac_manage_2pt_stripped.cc \
    formfac/formfac_manage_E_key.cc \
    formfac/formfac_manage_E_stripped.cc \
    formfac/formfac_manage_Z_key.cc \
    formfac/formfac_manage_Z_stripped.cc \
    formfac/formfac_manage_3pt.cc \
    formfac/formfac_manage_3pt_key.cc \
    formfac/formfac_manage_3pt_db.cc \
    formfac/formfac_manage_aggregate.cc \
    formfac/formfac_data_covar.cc \
    formfac/formfac_ensem_linfit.cc \
    formfac/formfac_linear_covar_genfit.cc \
    formfac/formfac_chisqq.cc \
    formfac/hadron_1pt_corr.cc \
    formfac/hadron_2pt_corr.cc \
    formfac/hadron_3pt_corr.cc

# Parton distribution sources
libadat_a_SOURCES += \
    parton/parton_correlation_functions_1.cc \
    parton/parton_correlation_functions_2.cc \
    parton/parton_derivatives.cc parton/parton_distribution_moments.cc

