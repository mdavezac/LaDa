include ../Files
bin_PROGRAMS = lamarck
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
	      @MPI_CFLAGS@ @BOOST_CXXFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@
lamarck_SOURCES = main.cc 
lamarck_CPPFLAGS =  $(AM_CPPFLAGS) -I .. 
lamarck_LDFLAGS =  $(AM_LDFLAGS)  @MPI_LIBS@ @BOOST_LDFLAGS@ @BOOST_MPI@ @BOOST_PROGOP@
lamarck_LDADD = liblamarck.la \
		../atat/libatat.la \
		../physics/libphysics.la \
		@MPI_ADDLIB@ \
		../print/libprint.la 
noinst_LTLIBRARIES = liblamarck.la
liblamarck_la_SOURCES = $(CRYSTAL_CC) $(CE_CC)
liblamarck_ladir = $(includedir)/lamarck
nobase_liblamarck_la_HEADERS = $(CRYSTAL_H) $(CE_H)
liblamarck_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
../atat/libatat.la: 
	cd ../atat; $(MAKE) 

../atat/libphysic.la: 
	cd ../physics; $(MAKE) 

../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)

../print/libprint.la:
	cd ../print; $(MAKE)
