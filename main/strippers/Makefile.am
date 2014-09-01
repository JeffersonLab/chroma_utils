#
# Boiler Plate stuff: Autoconf Substitutions etc
#
AM_CXXFLAGS = -I@top_srcdir@/lib -I@top_builddir@/include \
              -I@top_srcdir@/other_libs/xpath_reader/include \
              -I@top_srcdir@/other_libs/filedb/src \
              -I@top_srcdir@/other_libs/filedb/filehash \
              @CXXFLAGS@ @LIBXML2_CXXFLAGS@

AM_LDFLAGS = -L@top_builddir@/lib \
             -L@top_builddir@/other_libs/xpath_reader/lib \
             -L@top_builddir@/other_libs/filedb/src \
             -L@top_builddir@/other_libs/filedb/filehash \
             @LDFLAGS@

LDADD = -ladat -lXPathReader -lxmlWriter -lfiledb -lfilehash @LIBXML2_LIBS@ @LIBS@

#
# Local Headers
#
HDRS = strippers.h

#
# The programs to build
# 
## NOTE: if you want the programs to be build by default when 'make' is run
##   add them to bin_PROGRAMS.  But, if you don't want them to be installed
##   in @prefix@/bin by 'make install' then add them to nodist_PROGRAMS.
##   Another option is to add it to check_PROGRAMS, which means it will only
##   be built by 'make check' and is not installed.
bin_PROGRAMS = strip_spectrum strip_bar3ptfn strip_mresZa strip_dwfloc \
      strip_wallformfac strip_mres strip_mres4D strip_buildingblocks \
      strip_spectrumOct strip_spectrumQll strip_nCG get_pions strip_mesonspec \
      strip_hadspec strip_staticlight

check_PROGRAMS = 

EXTRA_PROGRAMS = strip_npr_vertex strip_hadron_contract

#
# The program and its dependencies
#
strip_hadspec_SOURCES =  strip_hadspec.cc $(HDRS)
strip_staticlight_SOURCES =  strip_staticlight.cc $(HDRS)
strip_mesonspec_SOURCES =  strip_mesonspec.cc $(HDRS)
strip_spectrum_SOURCES =  strip_spectrum.cc $(HDRS)
strip_spectrumOct_SOURCES =  strip_spectrumOct.cc $(HDRS)
strip_spectrumQll_SOURCES =  strip_spectrumQll.cc $(HDRS)
strip_bar3ptfn_SOURCES =  strip_bar3ptfn.cc $(HDRS)
strip_mresZa_SOURCES =  strip_mresZa.cc $(HDRS)
strip_mres_SOURCES =  strip_mres.cc $(HDRS)
strip_mres4D_SOURCES =  strip_mres4D.cc $(HDRS)
strip_nCG_SOURCES = strip_nCG.cc $(HDRS)
strip_dwfloc_SOURCES =  strip_dwfloc.cc $(HDRS)
strip_wallformfac_SOURCES =  strip_wallformfac.cc $(HDRS)
strip_buildingblocks_SOURCES =  strip_buildingblocks.cc $(HDRS)
strip_npr_vertex_SOURCES =  strip_npr_vertex.cc $(HDRS)
get_pions_SOURCES = get_pions.cc $(HDRS)
strip_hadron_contract_SOURCES =  strip_hadron_contract.cc hadron_contract_aggregate.cc \
   hadron_contract.cc hadron_2pt.cc simple_meson_2pt_w.cc 

#
# The folow rule will always try to rebuild libadat.a when you 
# try to compile spectrum
#
DEPENDENCIES = build_lib
#${noinst_PROGRAMS}: ${DEPENDENCIES}
${bin_PROGRAMS}: ${DEPENDENCIES}
##${check_PROGRAMS}: ${DEPENDENCIES}  # uncomment this when there are checkprograms
##${EXTRA_PROGRAMS}: ${DEPENDENCIES}

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_lib
build_lib:
	cd @top_builddir@/lib ; $(MAKE)
