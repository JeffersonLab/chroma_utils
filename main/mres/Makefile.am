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
#VPATH = @top_srcdir@/main/mres
#VPATH += @top_srcdir@/include

#
# Local Headers
#
HDRS=mres.h 

#
# The programs to build
# 
## NOTE: if you want the programs to be build by default when 'make' is run
##   add them to bin_PROGRAMS.  But, if you don't want them to be installed
##   in @prefix@/bin by 'make install' then add them to nodist_PROGRAMS.
##   Another option is to add it to check_PROGRAMS, which means it will only
##   be built by 'make check' and is not installed.
#check_PROGRAMS = 
#noinst_PROGRAMS = mres calc_Za
bin_PROGRAMS = mres calc_Za

#
# The program and its dependencies
#
mres_SOURCES =  mres.cc $(HDRS)
mres_DEPENDENCIES = build_lib

calc_Za_SOURCES =  calc_Za.cc $(HDRS)
calc_Za_DEPENDENCIES = build_lib


#
# The folow rule will always try to rebuild libadat.a when you 
# try to compile PROGRAMS
#
DEPENDENCIES = build_lib
${bin_PROGRAMS}: ${DEPENDENCIES}   # uncomment this when there are binprograms
#${check_PROGRAMS}: ${DEPENDENCIES}
#${EXTRA_PROGRAMS}: ${DEPENDENCIES}

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_lib
build_lib:
	cd @top_builddir@/lib ; $(MAKE)
