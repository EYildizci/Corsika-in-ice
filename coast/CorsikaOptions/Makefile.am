###########################################################################
# Makefile template
# $Id: Makefile.am,v 1.5 2007-10-22 13:46:01 pierog Exp $
###########################################################################

srcdir = @srcdir@
VPATH = @srcdir@

AM_CXXFLAGS = 
nobase_dist_include_HEADERS = 
EXTRA_DIST =
SUBDIRS =

AM_CXXFLAGS += -O3

if COAST_INTERNAL 

if USE_ROOTOUT
  SUBDIRS +=  rootout 
endif

if USE_INCLINED
  SUBDIRS += InclinedPlane
endif

if USE_COREAS
  SUBDIRS += CoREAS
endif

DIST_SUBDIRS = rootout InclinedPlane CoREAS

#distclean-local:
#	rm -f *.a ; rm -f *.la ; rm -f *.so*

#dist-hook:
#	tar --remove-file -zcf bernlohr-1.46.tar.gz [^z]*.c *.h *1.dat *2.dat *3.dat *4.dat *5.dat *6.dat *7.dat *8.dat *9.dat *.ps *.pdf Copyright GNUmakefile.org README  ; 


# --------------------------------------   END OF INTERNAL COAST BUILD ------------------------------------------------------------------
else


endif
