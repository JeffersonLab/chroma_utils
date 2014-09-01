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
#VPATH = @top_srcdir@/main/utils
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
bin_PROGRAMS = merge_proplists merge_cproplists  average_proplist concat_and_scale_proplist concat_proplist fold shift_corrfunc jack_ratio time_reverse block_proplist eff_mass jack_corrfunc_ratio plot_prop average_cproplist cfold jack_fit jack_ratio2 chop_proplist


#
# The program and its dependencies
#
concat_proplist_SOURCES =  concat_proplist.cc
concat_proplist_CPPFLAGS =  
concat_proplist_DEPENDENCIES = build_lib

block_proplist_SOURCES =  block_proplist.cc
block_proplist_CPPFLAGS =  
block_proplist_DEPENDENCIES = build_lib

concat_and_scale_proplist_SOURCES =  concat_and_scale_proplist.cc
concat_and_scale_proplist_CPPFLAGS =  
concat_and_scale_proplist_DEPENDENCIES = build_lib

jack_fit_SOURCES =  jack_fit.cc
jack_fit_CPPFLAGS =  
jack_fit_DEPENDENCIES = build_lib

fold_SOURCES =  fold.cc
fold_CPPFLAGS =  
fold_DEPENDENCIES = build_lib

cfold_SOURCES =  cfold.cc
cfold_CPPFLAGS =  
cfold_DEPENDENCIES = build_lib

plot_prop_SOURCES =  plot_prop.cc
plot_prop_CPPFLAGS =  
plot_prop_DEPENDENCIES = build_lib

eff_mass_SOURCES =  eff_mass.cc
eff_mass_CPPFLAGS =  
eff_mass_DEPENDENCIES = build_lib

time_reverse_SOURCES =  time_reverse.cc
time_reverse_CPPFLAGS =  
time_reverse_DEPENDENCIES = build_lib

jack_ratio_SOURCES =  jack_ratio.cc
jack_ratio_CPPFLAGS =  
jack_ratio_DEPENDENCIES = build_lib

jack_ratio2_SOURCES =  jack_ratio2.cc
jack_ratio2_CPPFLAGS =  
jack_ratio2_DEPENDENCIES = build_lib

shift_corrfunc_SOURCES =  shift_corrfunc.cc
shift_corrfunc_CPPFLAGS =  
shift_corrfunc_DEPENDENCIES = build_lib

average_proplist_SOURCES =  average_proplist.cc
average_proplist_CPPFLAGS =  
average_proplist_DEPENDENCIES = build_lib

average_cproplist_SOURCES =  average_cproplist.cc
average_cproplist_CPPFLAGS =  
average_cproplist_DEPENDENCIES = build_lib

merge_proplists_SOURCES =  merge_proplists.cc
merge_proplists_CPPFLAGS =  
merge_proplists_DEPENDENCIES = build_lib

merge_cproplists_SOURCES =  merge_cproplists.cc
merge_cproplists_CPPFLAGS =  
merge_cproplists_DEPENDENCIES = build_lib

jack_corrfunc_ratio_SOURCES =  jack_corrfunc_ratio.cc
jack_corrfunc_ratio_CPPFLAGS =  
jack_corrfunc_ratio_DEPENDENCIES = build_lib

chop_proplist_SOURCES =  chop_proplist.cc
chop_proplist_CPPFLAGS =  
chop_proplist_DEPENDENCIES = build_lib

#
# The folow rule will always try to rebuild libadat.a when you 
# try to compile a main prog
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
