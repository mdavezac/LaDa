if DOPESCAN
  PDARWIN=pdarwin
endif
bin_PROGRAMS = cedarwin $(PDARWIN)
AM_CPPFLAGS = -malign-double @ATAT_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ @MPICH_CFLAGS@ @EO_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@ @EO_LIBS@
cedarwin_SOURCES = main.cc ce.cc
cedarwindir = darwin
nobase_cedarwin_HEADERS = ce.h
cedarwin_CPPFLAGS =  -D _CE $(AM_CPPFLAGS) -I .. @GSL_INCLUDE@ 
cedarwin_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIB@ @MPICH_LIBS@
cedarwin_LDADD = libdarwin.la ../lamarck/liblamarck.la ../atat/libatat.la @MPI_ADDLIB@
pdarwin_SOURCES = main.cc pescan.cc
pdarwindir = darwin
nobase_pdarwin_HEADERS = pescan.h
pdarwin_CPPFLAGS =  -D _PESCAN $(AM_CPPFLAGS) -I .. @GSL_INCLUDE@ 
pdarwin_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIB@ @MPICH_LIBS@
pdarwin_LDADD = libdarwin.la ../lamarck/liblamarck.la ../atat/libatat.la ../vff/libvff.la \
	       	../pescan_interface/libpescan.la ../physics/libphysics.la @MPI_ADDLIB@
noinst_LTLIBRARIES = libdarwin.la
libdarwin_ladir = darwin
libdarwin_la_SOURCES = concentration.cc print_xmgrace.cc
nobase_libdarwin_la_HEADERS = breeder.h  checkpoints.h darwin.h evaluator.h gencount.h loadsave.h \
                              operators.h  print_xmgrace.h  statistics.h concentration.h darwin.impl.h \
			      functors.h individual.h minimizergenop.h results.h taboos.h
libopt_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
../atat/libatat.la: 
	cd ../atat; $(MAKE) 
../lamarck/liblamarck.la: 
	cd ../lamarck; $(MAKE)
../physics/libphysics.la:
	cd ../physics; $(MAKE)
../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)
../mpi/libpescan.la:
	cd ../pescan; $(MAKE)
