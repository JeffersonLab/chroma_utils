#
# Boiler Plate stuff: Autoconf Substitutions etc
#
AM_CXXFLAGS = -I@top_srcdir@/lib \
              -I@top_srcdir@/other_libs/xpath_reader/include \
              -I@top_srcdir@/other_libs/filedb/src \
              -I@top_srcdir@/other_libs/filedb/filehash \
              -I@top_srcdir@/include/recipes \
              @CXXFLAGS@ @LIBXML2_CXXFLAGS@

AM_LDFLAGS = -L. \
             -L@top_builddir@/lib \
             -L@top_builddir@/other_libs/xpath_reader/lib \
             -L@top_builddir@/other_libs/filedb/src \
             -L@top_builddir@/other_libs/filedb/filehash \
             @LDFLAGS@

LDADD = -lformfac -ladat -lXPathReader -lxmlWriter -lfiledb -lfilehash @LIBXML2_LIBS@ @LIBS@

#
# Library
#
lib_LIBRARIES = libformfac.a

libformfac_a_SOURCES = \
    matrix_elems.cc meson_matrix_elems.cc wavefuncs.cc operators.cc \
    current_ops.cc 


#
# The programs to build
# 
## NOTE: if you want the programs to be build by default when 'make' is run
##   add them to bin_PROGRAMS.  But, if you don't want them to be installed
##   in @prefix@/bin by 'make install' then add them to nodist_PROGRAMS.
##   Another option is to add it to check_PROGRAMS, which means it will only
##   be built by 'make check' and is not installed.

#noinst_PROGRAMS = spectrum

bin_PROGRAMS =  meson_3pt_fit_proj_snk thin_bb build_cache 
noinst_PROGRAMS = file_test make_mom_list

#check_PROGRAMS = 
EXTRA_PROGRAMS = nucleon_ff eigen_inverse complex_eigen_inverse test_manage


#
# The program and its dependencies
#
meson_3pt_fit_proj_snk_SOURCES = meson_3pt_fit_proj_snk.cc

eigen_inverse_SOURCES = eigen_inverse.cc
complex_eigen_inverse_SOURCES = complex_eigen_inverse.cc

thin_bb_SOURCES =  thin_bb.cc
build_cache_SOURCES =  build_cache.cc

#nucleon_ff_SOURCES =  nucleon_ff.cc ff_driver.cc state.h

#pion_src_pion_proj_snk_SOURCES = pion_src_pion_proj_snk.cc ff_driver.cc state.h
#rho_src_pion_proj_snk_SOURCES = rho_src_pion_proj_snk.cc ff_driver.cc state.h

test_manage_SOURCES = test_manage.cc

file_test_SOURCES =  file_test.cc

make_mom_list_SOURCES =  make_mom_list.cc

#
# The folow rule will always try to rebuild libadat.a when you 
# try to compile spectrum
#
DEPENDENCIES = build_lib ${lib_LIBRARIES}
${bin_PROGRAMS}: ${DEPENDENCIES}   # uncomment this when there are binprograms
${noinst_PROGRAMS}: ${DEPENDENCIES}
${check_PROGRAMS}: ${DEPENDENCIES}
${EXTRA_PROGRAMS}: ${DEPENDENCIES}

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_lib
build_lib:
	cd @top_builddir@/lib ; $(MAKE)
