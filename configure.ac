#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(corsika,74005)
AC_REVISION($Id: $)

AC_PREREQ(2.65)

########################################
# version 
offic="offic"     # official
#offic="noffic"     # not official


AC_CONFIG_HEADERS(include/config.h)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(src/corsika.F)
AC_PREFIX_DEFAULT(.)             #set default_prefix to . 
AC_PREFIX_PROGRAM()              #set prefix=default_prefix

AM_MAINTAINER_MODE
AC_CANONICAL_HOST

AM_INIT_AUTOMAKE([foreign subdir-objects])

dnl force compiler warnings
if ! echo $CXXFLAGS | grep -q -- -Wall ; then
  CXXFLAGS="$CXXFLAGS -Wall"
fi
if test "$FFLAGS" = "" ; then
  FFLAGS="-g"
fi
if test "$CFLAGS" = "" ; then
  CFLAGS="-g"
fi

# Checks for programs.
AC_PROG_F77(pgf77 ifc ifort gfortran f77 xlf fort77 g77)
AC_PROG_CC(cc gcc)
AC_PROG_CXX
AC_CHECK_PROG(CPP,cpp,"cpp")
AC_PROG_CPP
AC_PROG_LIBTOOL
AC_F77_LIBRARY_LDFLAGS
AC_CACHE_SAVE 
AM_PROG_CC_C_O

AX_BOOST_BASE
AX_BOOST_PROGRAM_OPTIONS   
AM_CONDITIONAL(HAVE_BOOST_PROGRAM_OPTIONS, test "x$ax_cv_boost_program_options" = "xyes")

CORSIKA_INIT

CORSIKA_ENABLE_DIST  #should be define first
CORSIKA_CONFIGURE_COMPILER
CORSIKA_ENABLE_M32
CORSIKA_ENABLE_OUT32 

CORSIKA_CONFIGURE_DETECTOR
CORSIKA_CONFIGURE_TIME
CORSIKA_CONFIGURE_HEMODEL
CORSIKA_CONFIGURE_LEMODEL

#thin should be before for proper name definition
CORSIKA_ENABLE_THIN
CORSIKA_ENABLE_CERENKOV
CORSIKA_ENABLE_IACT
CORSIKA_ENABLE_INTCLONG
CORSIKA_ENABLE_NOCLONG
CORSIKA_ENABLE_CERWLEN
CORSIKA_ENABLE_CEFFIC
CORSIKA_ENABLE_IACTEXT
CORSIKA_ENABLE_ATMEXT
CORSIKA_ENABLE_NEUTRINO
CORSIKA_ENABLE_PLOTSH
CORSIKA_ENABLE_PLOTSH2
CORSIKA_ENABLE_INTTEST
CORSIKA_ENABLE_SLANT
CORSIKA_ENABLE_CURVED
CORSIKA_ENABLE_UPWARD
CORSIKA_ENABLE_VIEWCONE
CORSIKA_ENABLE_AUGCERLONG
CORSIKA_ENABLE_PRESHOWER
CORSIKA_ENABLE_MUPROD
CORSIKA_ENABLE_COMPACT
CORSIKA_ENABLE_ANNITEST
CORSIKA_ENABLE_LPM
CORSIKA_ENABLE_STACKIN
CORSIKA_ENABLE_CHARM
CORSIKA_ENABLE_HERWIG

CORSIKA_ENABLE_CHARM
CORSIKA_ENABLE_TAULEP
CORSIKA_ENABLE_EHISTORY
CORSIKA_ENABLE_TRAJECT
CORSIKA_ENABLE_EFIELD

CORSIKA_ENABLE_ROOTOUT
CORSIKA_ENABLE_MACHINEINDEPENDENT
CORSIKA_ENABLE_COASTUSERLIB
CORSIKA_ENABLE_COREAS
CORSIKA_ENABLE_INCLINED

#conex should be after coast for proper compilation 
CORSIKA_ENABLE_CONEX

#hbook output files should be after CONEX for proper name definitions
CORSIKA_ENABLE_ANAHIST
CORSIKA_ENABLE_AUGERINFO
CORSIKA_ENABLE_AUGERHIST
CORSIKA_ENABLE_MUONHIST

CORSIKA_ENABLE_PARALLEL
#at the end to be sure all libraries are included properly
CORSIKA_ENABLE_PARALLELIB
CORSIKA_ENABLE_NOCOMPILE
CORSIKA_ENABLE_KEEPSOURCE

CORSIKA_NAME="$PACKAGE$VERSION$USEROS$HIMODEL$LOMODEL$INTT"
CORSIKA_NAME="$PLIB1$CORSIKA_NAME$PLIB2"

AC_SUBST([BERNLIBS]) # list of bernlohr libraries to be compiled 

CORSIKA_LIBS="$CORSIKA_LIBS $COAST_INTERFACE_LIBS"

#
# These conditionals are needed to set the conditions to false when they are not needed
#
AM_CONDITIONAL([COAST_INTERNAL], test "true" = "true" ) # immer true
AM_CONDITIONAL([HAVE_DOXYGEN], test "x$COAST_DIR" = "xNODOXYGEN" ) # CORSIKA never with doxygen !!!
AM_CONDITIONAL([HAVE_ROOT], test x$corsika_cv_root_root != xno)

AC_SUBST([CORSIKA_NAME])
AC_SUBST([CORSIKA_SUBDIRS])
AC_SUBST([CORSIKA_LIB_DIR])
AC_SUBST([CORSIKA_LIBS])
AC_SUBST([MODEL_LIBS])
AC_SUBST([CORSIKA_FFLAGS])
AC_SUBST([MPIRUNNER_LIBS])

if test  "x$offic" = "xoffic" ; then
  AC_DEFINE([__OFFIC__], 1, [official version marker])
else
  AC_DEFINE([__OFFIC__], 0, [official version marker])
fi


if test  "x$corsika_cv_enable_DIST" != "xyes" ; then
AC_SUBST([HIMODEL])

CORSIKA_CHECK_MISSING_REQUIRED_PACKAGES
CORSIKA_CHECK_CONFLICTS
CORSIKA_CHECK_REQUIRED
fi


AC_CONFIG_FILES([Makefile  bernlohr/Makefile conex/Makefile
		dpmjet/Makefile epos/Makefile pythia/Makefile
                herwig/Makefile nexus/Makefile urqmd/Makefile src/Makefile 
		run/Makefile doc/Makefile lib/Makefile 
		coast/Makefile
		coast/Documentation/Makefile
		coast/CorsikaOptions/rootout/Makefile
		coast/CorsikaOptions/CoREAS/Makefile
		coast/CorsikaOptions/Makefile
		coast/CorsikaOptions/InclinedPlane/Makefile
		coast/CorsikaFileIO/Makefile
		coast/CorsikaInterface/Makefile
		coast/CorsikaToROOT/Makefile
		coast/CorsikaROOT/Makefile
		coast/CorsikaIntern/Makefile])


AC_OUTPUT

if test  "x$corsika_cv_enable_DIST" != "xyes" ; then
CORSIKA_FINISH
fi

