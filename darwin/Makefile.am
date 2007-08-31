if DOPESCAN
  PDARWIN=pdarwin
endif
if DOCE
  CEDARWIN= cedarwin
endif
bin_PROGRAMS = $(CEDARWIN) $(PDARWIN)
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ @MPICH_CFLAGS@ @EO_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@ @EO_LIBS@
cedarwin_SOURCES = main.cc ce.cc single_site.cc
cedarwindir = darwin
nobase_cedarwin_HEADERS = ce.h single_site.impl.h single_site.h
cedarwin_CPPFLAGS =  -D _CE $(AM_CPPFLAGS) -I .. @GSL_INCLUDE@ 
cedarwin_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIB@ @MPICH_LIBS@
cedarwin_LDADD = libdarwin.la ../lamarck/liblamarck.la ../atat/libatat.la\
		 ../print/libprint.la @MPI_ADDLIB@
pdarwin_SOURCES = main.cc pescan.cc two_sites.cc single_site.cc
pdarwin_HEADERS = pescan.h two_sites.h two_sites.impl.h single_site.h
pdarwindir = darwin
nobase_pdarwin_HEADERS = pescan.h
pdarwin_CPPFLAGS =  -D _PESCAN $(AM_CPPFLAGS) -I .. @GSL_INCLUDE@ 
pdarwin_LDFLAGS =  $(AM_LDFLAGS) @GSL_LIB@ @MPICH_LIBS@
pdarwin_LDADD = ../lamarck/liblamarck.la libdarwin.la ../lamarck/liblamarck.la ../atat/libatat.la ../vff/libvff.la \
	       	../pescan_interface/libpescan.la ../physics/libphysics.la ../print/libprint.la @MPI_ADDLIB@
noinst_LTLIBRARIES = libdarwin.la
libdarwin_ladir = darwin
libdarwin_la_SOURCES = concentration.cc objective.cc 
nobase_libdarwin_la_HEADERS = breeder.h concentration.h evaluation.h functors.h loadsave.h operators.h \
	                      statistics.h darwin.h gencount.h store.h \
			      checkpoints.h darwin.impl.h evaluator.h individual.h objective.h \
			      taboos.h gatraits.h 
libopt_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
../atat/libatat.la: 
	cd ../atat; $(MAKE) 
../lamarck/liblamarck.la: 
../physics/libphysics.la:
	cd ../physics; $(MAKE)
../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)
../pescan_interface/libpescan.la:
	cd ../pescan; $(MAKE)
../print/libprint.la:
	cd ../print; $(MAKE)
