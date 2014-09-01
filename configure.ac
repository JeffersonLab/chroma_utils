AC_INIT(adat,3.0,edwards@jlab.org)

AC_CONFIG_AUX_DIR(config)

AC_LANG_CPLUSPLUS
AC_PROG_CXX(g++ cxx CC)
AC_PROG_CC(gcc cc)
AC_PROG_RANLIB

dnl Initialise automake
AM_INIT_AUTOMAKE
dnl
dnl
dnl --with-xxxx and --enable-xxxx switches 
dnl

dnl
dnl libxml2
dnl
AC_ARG_WITH(libxml2,
  AC_HELP_STRING(
    [--with-libxml2=DIR],
    [Build on top of libxml2 where libxml2 is installed in DIR]
  ),
  [LIBXML2_HOME="$with_libxml2"]
)


echo ${LIBXML2_HOME}

if test "X${LIBXML2_HOME}X" = "XX" ; then
  AC_PATH_PROG(LIBXML2_CONFIG, [xml2-config], [])
else
  AC_PATH_PROG(LIBXML2_CONFIG, [xml2-config], [], [${LIBXML2_HOME}/bin:${PATH}])
fi

if test "X${LIBXML2_CONFIG}X" = "XX" ; then
  AC_MSG_ERROR([libxml2 configuration program xml2-config not found.])
fi

AC_MSG_NOTICE([Found libxml2 configuration program ${LIBXML2XX_CONFIG}])
AC_SUBST(LIBXML2_CXXFLAGS, "`${LIBXML2_CONFIG} --cflags`")
AC_MSG_NOTICE([libxml2 compile flags: ${LIBXML2_CXXFLAGS}])
AC_SUBST(LIBXML2_LIBS,     "`${LIBXML2_CONFIG} --libs`")
AC_MSG_NOTICE([libxml2 libraries flags: ${LIBXML2_LIBS}])

dnl This function is defined in acinclude.m4
dnl it tries to try to compile the program 
AC_MSG_CHECKING([if we can compile/link a simple libxml2 program])
PAC_LIBXML2_LINK_CXX_FUNC(
  ${LIBXML2_CXXFLAGS},
  ${LIBXML2_LIBS},
  ,
  ,
  [libxml2_link_ok=yes],
  [libxml2_link_ok=no]
)

if test "X${libxml2_link_ok}X" = "XyesX" ; then 
  AC_MSG_RESULT(yes)
else 
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([Cannot compile/link a program with libxml2.
    Use --with-libxml2=<dir> to select a working version.])
fi 


dnl itpp support
AC_ARG_WITH(itpp, 
  AC_HELP_STRING(
    [--with-itpp=DIR],
    [Use ITPP in DIR]),
    [ITPP_HOME="$with_itpp"],
    [ITPP_HOME="no"]
)


if test "X${ITPP_HOME}X" != "XnoX"; then
  AC_MSG_NOTICE("Using ITPP in ${ITPP_HOME}")
  AC_PATH_PROG(ITPP_CONFIG, [itpp-config], [], [${ITPP_HOME}/bin:$PATH] )
  if test "X${ITPP_CONFIG}X" = "XX"; then 
    AC_MSG_ERROR([Program itpp-config not found])
  fi
  AC_MSG_NOTICE([Found ITPP configuration program ${LIBXML2XX_CONFIG}])
  AC_SUBST(ITPP_CXXFLAGS, "`${ITPP_CONFIG} --cflags`")
  AC_MSG_NOTICE([ITPP compile flags: ${ITPP_CXXFLAGS}])
  AC_SUBST(ITPP_LIBS,     "`${ITPP_CONFIG} --libs`")
  AC_MSG_NOTICE([ITPP libraries flags: ${ITPP_LIBS}])
fi


dnl Enable mpi support
AC_ARG_ENABLE(mpi,
   AC_HELP_STRING([--enable-mpi], 
     [Enable MPI support in some routines]),
     [ac_mpi_enabled=1
      AC_MSG_NOTICE([Enabling MPI support in some routines])],
     [ac_mpi_enabled=0]
)

if test ${ac_mpi_enabled} -eq 1; then
   AC_DEFINE_UNQUOTED(USE_MPI, ${ac_mpi_enabled}, [Enable MPI support])
fi

AM_CONDITIONAL(ADAT_MPI_COMMS,    [test "x${ac_mpi_enabled}x" = "x1x" ])
AM_CONDITIONAL(ADAT_SINGLE_COMMS, [test "x${ac_mpi_enabled}x" = "x0x" ])


dnl
dnl
dnl Now have all the options... Do some configuration 
dnl
dnl
##################################
# Check for programs to build docs
##################################
AC_CHECK_PROG(LATEX, latex, latex, true)
AC_CHECK_PROG(DVIPS, dvips, dvips, true)
AC_CHECK_PROG(DVIPDF, dvipdf, dvipdf, true)
AC_CHECK_PROG(DOXYGEN, doxygen, doxygen, true)

#######################
# Produce output
#######################

#
# Configure subdir packages
#
AC_CONFIG_SUBDIRS(other_libs/xpath_reader)
AC_CONFIG_SUBDIRS(other_libs/minuit)
AC_CONFIG_SUBDIRS(other_libs/filedb)

#
# Set up the doxygen file
# Substitute the right directories
#
# AC_CONFIG_FILES(docs/adatdoc)

#
# Generate makefiles
#
AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(adat-config)
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(other_libs/Makefile)
# AC_CONFIG_FILES(docs/Makefile)
AC_CONFIG_FILES(main/Makefile)
AC_CONFIG_FILES(main/strippers/Makefile)
#AC_CONFIG_FILES(main/spectrum/Makefile)
#AC_CONFIG_FILES(main/mres/Makefile)
#AC_CONFIG_FILES(main/fpi/Makefile)
AC_CONFIG_FILES(main/fitters/Makefile)
AC_CONFIG_FILES(main/utils/Makefile)
#AC_CONFIG_FILES(main/formfac/Makefile)
#
# Finish
#
AC_OUTPUT

