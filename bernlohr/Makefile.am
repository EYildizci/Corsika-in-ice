## Makefile.am -- Process this file with automake to get Makefile.in

# Maximum output I/O buffer size, limiting the number of photon
# bunches per telescope that can be handled.
srcdir = @srcdir@
VPATH = @srcdir@

MAX_IO_BUFFER=200000000


EXTRA_LIBRARIES = libbern.a libiact.a libtestlib.a libatmext.a ## ?? needed ??

lib_LIBRARIES = @BERNLIBS@
libbern_a_SOURCES = fileopen.c straux.c
libbern_a_CFLAGS = @BERN_CFLAGS@

libiact_a_SOURCES = iact.c io_simtel.c eventio.c warning.c sampling.c
libiact_a_CPPFLAGS = @BERNNOATMEXT_CPPFLAGS@ @BERNIACTEXT_CPPFLAGS@ -DMAX_IO_BUFFER=$(MAX_IO_BUFFER) \
		     -DCORSIKA_VERSION=$(PACKAGE_VERSION) 
libiact_a_CFLAGS = @BERN_CFLAGS@

libtestlib_a_SOURCES = ztestlib.c
libtestlib_a_CPPFLAGS = @BERN_CPPFLAGS@

libatmext_a_SOURCES = atmo.c fileopen.c straux.c
libatmext_a_CFLAGS = @BERN_CFLAGS@ 

EXTRA_DIST = Copyright GNUmakefile.org README README.ps atmo.h  \
	     atmprof1.dat atmprof2.dat atmprof3.dat atmprof4.dat \
	     atmprof5.dat atmprof6.dat atmprof7.dat atmprof8.dat \
             atmprof9.dat corsika_build_script gen_config trapfpe.c\
	     eventio_en.pdf fileopen.h iact3d.ps gpl.txt lgpl-2.1.txt lgpl-3.txt\
	     iact_refman.pdf initial.h io_basic.h mc_tel.h sampling.h \
	     sim_skeleton.c straux.h testio.c warning.h listio.c \
	     iact_refman.pdf ztestlib.c README.corsika

distclean-local:
	rm -f *.a ; rm -f *.la ; rm -f *.so*

README:
	@if test -f bernlohr-1.48.tar.gz ; then \
	tar zxf bernlohr-1.48.tar.gz ; rm -f bernlohr-1.48.tar.gz ;\
	fi

dist-hook:
	tar --remove-file -zcf bernlohr-1.48.tar.gz [^z]*.c *.h *1.dat *2.dat *3.dat *4.dat *5.dat *6.dat *7.dat *8.dat *9.dat *.ps *.pdf *.txt gen_config corsika_build_script Copyright GNUmakefile.org README  ; 
