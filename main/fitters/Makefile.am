#
# Boiler Plate stuff: Autoconf Substitutions etc
#
INCFLAGS = -I@top_srcdir@/lib -I@top_builddir@/include \
           -I@top_srcdir@/other_libs/xpath_reader/include

AM_CXXFLAGS = $(INCFLAGS) @CXXFLAGS@ @LIBXML2_CXXFLAGS@
AM_LDFLAGS = -L@top_builddir@/lib \
             -L@top_builddir@/other_libs/xpath_reader/lib \
             @LDFLAGS@
LDADD = -ladat -lXPathReader -lxmlWriter @LIBXML2_LIBS@ @LIBS@

#
# VPATH
#
#VPATH = @top_srcdir@/main/fitters
#VPATH += @top_srcdir@/include

#
# Local Headers
#
HDRS=

#
# The programs to build
# 
## NOTE: if you want the programs to be build by default when 'make' is run
##   add them to bin_PROGRAMS.  But, if you don't want them to be installed
##   in @prefix@/bin by 'make install' then add them to nodist_PROGRAMS.
##   Another option is to add it to check_PROGRAMS, which means it will only
##   be built by 'make check' and is not installed.
#check_PROGRAMS = 
#noinst_PROGRAMS = polyfit covarprop readcovp exponfit exponperfit cov_polyfit cov_exponfit cov_chiPTfpi chiPTfpi chiPTlog cov_chiPTlog alt_exponfit cov_alt_exponfit cov_exponperfit
bin_PROGRAMS = polyfit covarprop readcovp exponfit exponperfit cov_polyfit cov_exponfit cov_chiPTfpi chiPTfpi chiPTlog cov_chiPTlog alt_exponfit cov_alt_exponfit cov_exponperfit boundmom cov_boundmom
#
# The program and its dependencies
#
polyfit_SOURCES =  fit.cc ../../include/covfit/fitters/fitfunctions.h
polyfit_CPPFLAGS =  -DFUNCTION=Polynomium
polyfit_DEPENDENCIES = build_lib

boundmom_SOURCES =  fit.cc ../../include/covfit/fitters/fitfunctions.h
boundmom_CPPFLAGS =  -DFUNCTION=BoundMom
boundmom_DEPENDENCIES = build_lib

exponfit_SOURCES =  fit.cc ../../include/covfit/fitters/fitfunctions.h 
exponfit_CPPFLAGS =  -DFUNCTION=Exponentials
exponfit_DEPENDENCIES = build_lib


alt_exponfit_SOURCES =  fit.cc ../../include/covfit/fitters/fitfunctions.h 
alt_exponfit_CPPFLAGS =  -DFUNCTION=AltExponentials
alt_exponfit_DEPENDENCIES = build_lib

exponperfit_SOURCES =  fit.cc ../../include/covfit/fitters/fitfunctions.h
exponperfit_CPPFLAGS =  -DFUNCTION=ExponentialsPeriodic
exponperfit_DEPENDENCIES = build_lib

cov_polyfit_SOURCES =  cofit.cc ../../include/covfit/fitters/fitfunctions.h
cov_polyfit_CPPFLAGS =  -DFUNCTION=Polynomium
cov_polyfit_DEPENDENCIES = build_lib

cov_boundmom_SOURCES =  cofit.cc ../../include/covfit/fitters/fitfunctions.h
cov_boundmom_CPPFLAGS =  -DFUNCTION=BoundMom
cov_boundmom_DEPENDENCIES = build_lib

cov_exponfit_SOURCES =  cofit.cc ../../include/covfit/fitters/fitfunctions.h 
cov_exponfit_CPPFLAGS =  -DFUNCTION=Exponentials
cov_exponfit_DEPENDENCIES = build_lib

cov_exponperfit_SOURCES =  cofit.cc ../../include/covfit/fitters/fitfunctions.h
cov_exponperfit_CPPFLAGS =  -DFUNCTION=ExponentialsPeriodic
cov_exponperfit_DEPENDENCIES = build_lib


cov_alt_exponfit_SOURCES =  cofit.cc ../../include/covfit/fitters/fitfunctions.h 
cov_alt_exponfit_CPPFLAGS =  -DFUNCTION=AltExponentials
cov_alt_exponfit_DEPENDENCIES = build_lib

chiPTfpi_SOURCES =  fit.cc ../../include/covfit/fitters/fitfunctions.h 
chiPTfpi_CPPFLAGS =  -DFUNCTION=ChiPTfpi
chiPTfpi_DEPENDENCIES = build_lib

cov_chiPTfpi_SOURCES =  cofit.cc ../../include/covfit/fitters/fitfunctions.h 
cov_chiPTfpi_CPPFLAGS =  -DFUNCTION=ChiPTfpi
cov_chiPTfpi_DEPENDENCIES = build_lib

chiPTlog_SOURCES =  fit.cc ../../include/covfit/fitters/fitfunctions.h 
chiPTlog_CPPFLAGS =  -DFUNCTION=ChiPTlog
chiPTlog_DEPENDENCIES = build_lib

cov_chiPTlog_SOURCES =  cofit.cc ../../include/covfit/fitters/fitfunctions.h 
cov_chiPTlog_CPPFLAGS =  -DFUNCTION=ChiPTlog
cov_chiPTlog_DEPENDENCIES = build_lib

#
# The latter rule will always try to rebuild libadat.a when you 
# try to compile a polyfit

covarprop_SOURCES =  covarprop.cc 
covarprop_DEPENDENCIES = build_lib

readcovp_SOURCES =  readcovp.cc 
readcovp_DEPENDENCIES = build_lib


# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_lib
build_lib:
	cd @top_builddir@/lib ; $(MAKE)
