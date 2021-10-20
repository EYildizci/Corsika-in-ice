###########################################################################
# Makefile template
# $Id: Makefile.am,v 1.2 2007-10-19 11:13:40 rulrich Exp $
###########################################################################

srcdir = @srcdir@
VPATH = @srcdir@

AM_CXXFLAGS = 
nobase_dist_include_HEADERS = 
EXTRA_DIST =



# ACLOCAL_AMFLAGS = -I m4
# 
#nobase_nodist_include_HEADERS = 

# documentation
#docdir = $(pkgdatadir)/doc
#nobase_doc_DATA = 
#nobase_nodist_doc_DATA =
#nobase_dist_doc_DATA =




if COAST_INTERNAL 

AM_CXXFLAGS += -I$(top_builddir)/coast/CorsikaIntern

# --------------------------------------   END OF INTERNAL COAST BUILD ------------------------------------------------------------------
else

AM_CXXFLAGS += -I$(top_builddir)/CorsikaIntern

# coast_lib_dir = $(top_builddir)

endif

###########################################################################
# globals
###########################################################################


AM_CXXFLAGS += -O3

###########################################################################
# crsRead
###########################################################################

crsRead_Headers =
noinst_crsRead_Headers =

crsRead:
	mkdir -p $@

crsRead-stamp: $(crsRead_Headers) $(noinst_crsRead_Headers) Makefile
	for i in $^; do \
	  ln -fs @abs_srcdir@/`basename $$i` crsRead/`basename $$i`; \
	done
	touch $@

crsRead-clean:
	rm -f crsRead-stamp
	rm -rf crsRead

crsRead_Headers += TCorsikaReader.h TSubBlockIO.h \
                   MParticleSubBlockOutput.h  \
                   MCorsikaReader.h  \
                   MLongfileReader.h

BUILT_SOURCES = crsRead crsRead-stamp


if COAST_INTERNAL  # to prevent a general boost dependency for CORSIKA

lib_LIBRARIES = libCorsikaFileIO.a
libCorsikaFileIO_a_SOURCES = TCorsikaReader.cc TSubBlockIO.cc \
                              MParticleSubBlockOutput.cc  \
                              MCorsikaReader.cc MLongfileReader.cc

libCorsikaFileIO_a_CPPFLAGS = #-fPIC
#libCorsikaFileIO_a_LDFLAGS =  --no-undefined -fPIC --shared

else

lib_LTLIBRARIES = libCorsikaFileIO.la
libCorsikaFileIO_la_SOURCES = TCorsikaReader.cc TSubBlockIO.cc \
                              MParticleSubBlockOutput.cc  \
                              MCorsikaReader.cc MLongfileReader.cc

libCorsikaFileIO_la_CPPFLAGS = -fPIC
libCorsikaFileIO_la_LDFLAGS = -fPIC --shared --no-undefined

bin_PROGRAMS = merge_corsika check_corsika
merge_corsika_SOURCES = merge_corsika.cc
merge_corsika_LDADD = @COAST_LIBS@
merge_corsika_LDFLAGS = -lboost_program_options 

check_corsika_SOURCES = check_corsika.cc
check_corsika_LDADD = @COAST_LIBS@
check_corsika_LDFLAGS = -lboost_program_options 

endif

#distclean-local:
#	rm -f *.a ; rm -f *.la ; rm -f *.so*; rm -rf crsRead-clean crsRead

clean-local: crsRead-clean

#install : $(BUILT_SOURCES) install-am

nobase_nodist_include_HEADERS = $(addprefix crsRead/, $(notdir $(crsRead_Headers)))
EXTRA_DIST += $(crsRead_Headers)
