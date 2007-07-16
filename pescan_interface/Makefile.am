AM_CPPFLAGS = -malign-double @ATAT_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ @MPICH_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@ 
noinst_LTLIBRARIES = libpescan.la
libpescan_la_SOURCES = interface.cc
libpescan_ladir = pescan_interface
nobase_libpescan_la_HEADERS = interface.h
libpescan_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
../atat/libatat.la: 
	cd ../atat; $(MAKE) 
../lamarck/liblamarck.la: 
	cd ../lamarck; $(MAKE)
../physics/libphysics.la:
	cd ../physics; $(MAKE)
../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)
