## Makefile.am -- Process this file with automake to get Makefile.in

SUBDIRS = @CORSIKA_SUBDIRS@

DIST_SUBDIRS = bernlohr coast conex dpmjet epos herwig pythia nexus urqmd doc run lib src .
EXTRA_DIST = coconut cDO_NOT_RUN_CONFIGURE
ACLOCAL_AMFLAGS = -I m4

if USE_PARALLELIB
install-exec-hook : 
	@cp -p -f @bindir@/mpi_runner$(EXEEXT) @bindir@/@CORSIKA_NAME@ ; \
	rm -f @bindir@/mpi_runner$(EXEEXT) ; \
	rm -f @prefix@/src/mpi_runner$(EXEEXT) ; \
	echo " " ;  \
	echo " " ;  \
	echo " " ;  \
	echo '--> "@CORSIKA_NAME@" successfully installed in :' ; \
	echo "    @bindir@/" ; \
	echo " " ; \
	echo " " 
else
install-exec-hook : 
	@cp -p -f @bindir@/corsika$(EXEEXT) @bindir@/@CORSIKA_NAME@ ; \
	rm -f @bindir@/corsika$(EXEEXT) ; \
	rm -f @prefix@/src/corsika$(EXEEXT) ; \
	echo " " ;  \
	echo " " ;  \
	echo " " ;  \
	echo '--> "@CORSIKA_NAME@" successfully installed in :' ; \
	echo "    @bindir@/" ; \
	echo " " ;  \
	if test "@HIMODEL@" = "_NEXUS" ; then \
	echo "--> You can run CORSIKA in @bindir@/ using for instance :" ; \
	echo "    ./@CORSIKA_NAME@ < all-inputs-nexus > output.txt" ; \
	elif test "@HIMODEL@" = "_EPOS" ; then \
	echo "--> You can run CORSIKA in @bindir@/ using for instance :" ; \
	echo "    ./@CORSIKA_NAME@ < all-inputs-epos > output.txt" ; \
	else \
	echo "--> You can run CORSIKA in @bindir@/ using for instance :" ; \
	echo "    ./@CORSIKA_NAME@ < all-inputs > output.txt" ; \
	fi ; \
	echo " " ;  \
	echo " " ;  \
	echo " " 
endif

if USE_NOCOMPILE

preprocess: 
	@cd src && $(MAKE) $(AM_MAKEFLAGS) nocompile ; \
	echo " " ;  \
	echo " " ;  \
	echo '--> To compile corsika, use "coconut -i".' ; \
	echo " " ;  \
	echo " " 

endif

dist-hook:
	cd $(distdir)/conex ; tar --remove-file -zcf conex-4.3700.tar.gz conex* egs4.* ===LOG=== README  ; cd ../..    ; \
	cd $(distdir)/epos ; tar --remove-file -zcf epos-lhc-v3400_cors.tar.gz epos*.f epos.inc* epos.ini* epos.param  Makefile.ka README ; cd ../..    ; \
	cd $(distdir)/nexus ; tar --remove-file -zcf nexus3.97_cors.tar.gz nexus*.f nexus.inc* nexus.ini* nexus.param Makefile.ka README ; cd ../..    ; \
	cd $(distdir)/urqmd ; tar --remove-file -zcf urqmd1.3_cors.tar.gz *.f Copyright GNUmakefile_corsika README ; cd ../..    ; \
	cd $(distdir)/bernlohr ; tar --remove-file -zcf bernlohr-1.48.tar.gz [^z]*.c *.h *1.dat *2.dat *3.dat *4.dat *5.dat *6.dat *7.dat *8.dat *9.dat *.ps *.pdf *.txt gen_config corsika_build_script Copyright GNUmakefile.org README  ; cd ../..    ; 
	cd $(distdir)/pythia ; tar --remove-file -zcf pythia-6.4.11.tar.gz README update* *.f *.pdf  ; cd ../.. 
#	rm -rf `find $(distdir)/doc -path '*svn*' ` ;
