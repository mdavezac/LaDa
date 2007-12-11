if DOBANDGAPOPT
  BANDGAPOPT=bandgap_opt
endif
if DOCEOPT
  CEOPT = ce_opt
endif
if DOLAYEREDOPT
  LAYEREDOPT = layered_opt
endif
if DOEMASSOPT
  EMASSOPT = emass_opt
endif
bin_PROGRAMS = $(CEOPT) $(BANDGAPOPT) $(LAYEREDOPT) $(EMASSOPT)
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ @MPICH_CFLAGS@ @EO_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@ @EO_LIBS@
ce_opt_SOURCES = main.cc ce.cc single_site.cc groundstate.cc
ce_opt_HEADERS = ce.h single_site.h single_site.impl.h groundstate.h
ce_optdir = $(includedir)/darwin
nobase_ce_opt_HEADERS = ce.h single_site.impl.h single_site.h groundstate.h
ce_opt_CPPFLAGS =  -D _CE $(AM_CPPFLAGS) -I .. @GSL_INCLUDE@ 
ce_opt_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIB@ @MPICH_LIBS@
ce_opt_LDADD = libdarwin.la ../lamarck/liblamarck.la ../physics/libphysics.la ../atat/libatat.la \
	       ../print/libprint.la @MPI_ADDLIB@ libdarwin.la 

bandgap_opt_SOURCES = main.cc bandgap_stubs.cc two_sites.cc single_site.cc vff.cc bandgap.cc 
bandgap_opt_HEADERS = bandgap_stubs.h two_sites.h two_sites.impl.h single_site.h \
		 single_site.impl.h vff.h bandgap.h 
bandgap_optdir = $(includedir)/darwin
bandgap_opt_CPPFLAGS =  -D _PESCAN $(AM_CPPFLAGS) -I .. @GSL_INCLUDE@ 
bandgap_opt_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIB@ @MPICH_LIBS@
bandgap_opt_LDADD = libdarwin.la ../lamarck/liblamarck.la ../atat/libatat.la ../vff/libvff.la \
	            ../pescan_interface/libpescan.la ../physics/libphysics.la ../print/libprint.la @MPI_ADDLIB@

layered_opt_SOURCES = main.cc bandgap_stubs.cc two_sites.cc single_site.cc vff.cc molecularity.cc layered.cc
layered_opt_HEADERS = bandgap_stubs.h bandgap_stubs.impl.h two_sites.h two_sites.impl.h \
		      single_site.h single_site.impl.h vff.h molecularity.h molecularity.impl.h \
		      layered.h layered.impl.h
layered_optdir = $(includedir)/darwin
layered_opt_CPPFLAGS =  -D _MOLECULARITY $(AM_CPPFLAGS) -I .. @GSL_INCLUDE@ 
layered_opt_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIB@ @MPICH_LIBS@
layered_opt_LDADD = libdarwin.la ../lamarck/liblamarck.la ../atat/libatat.la ../vff/libvff.la \
		    ../pescan_interface/libpescan.la ../physics/libphysics.la\
		    ../print/libprint.la @MPI_ADDLIB@

emass_opt_SOURCES = main.cc emass_stubs.cc two_sites.cc single_site.cc vff.cc emass.cc layered.cc
emass_opt_HEADERS = emass_stubs.h emass_stubs.impl.h two_sites.h two_sites.impl.h \
		    single_site.h single_site.impl.h vff.h emass.h emass.impl.h \
		    layered.h layered.impl.h
emass_optdir = $(includedir)/darwin
emass_opt_CPPFLAGS =  -D _EMASS $(AM_CPPFLAGS) -I .. @GSL_INCLUDE@ 
emass_opt_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIB@ @MPICH_LIBS@
emass_opt_LDADD = libdarwin.la ../lamarck/liblamarck.la ../atat/libatat.la ../vff/libvff.la \
		  ../pescan_interface/libpescan.la ../physics/libphysics.la \
		  ../print/libprint.la ../lapack/liblapack.la @MPI_ADDLIB@ @FCLIBS@

noinst_LTLIBRARIES = libdarwin.la
libdarwin_ladir = $(includedir)/darwin
libdarwin_la_SOURCES = concentration.cc debug.cc
nobase_libdarwin_la_HEADERS = breeder.h concentration.h evaluation.h functors.h loadsave.h \
			      operators.h statistics.h darwin.h gencount.h store.h \
			      checkpoints.h darwin.impl.h evaluator.h individual.h objective.h \
			      taboos.h gatraits.h  bitstring.h fitness.h gaoperators.h \
			      scaling.h darwin.impl.h  evaluation.impl.h individual.impl.h \
			      scaling.impl.h store.impl.h  debug.h fitness.impl.h \
			      taboos.impl.h objective.impl.h \
			      bitstring.impl.h gaoperators.impl.h checkpoints.impl.h \
			      minimizergenop.h minimizergenop.impl.h
libdarwin_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
../atat/libatat.la: 
	cd ../atat; $(MAKE) 
../lamarck/liblamarck.la: 
	cd ../lamarck; $(MAKE)
../physics/libphysics.la:
	cd ../physics; $(MAKE)
../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)
../pescan_interface/libpescan.la:
	cd ../pescan; $(MAKE)
../print/libprint.la:
	cd ../print; $(MAKE)
../print/liblapack.la:
	cd ../lapack; $(MAKE)
