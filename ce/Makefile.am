#
# Version: $Id$
#
include ../Files
bin_PROGRAMS = regular # cubic tetra 
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
	      @MPI_CFLAGS@ @BOOST_CXXFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@
cubic_SOURCES = main.cc 
cubic_CPPFLAGS =  $(AM_CPPFLAGS) -I .. -D_CUBIC_CE_
cubic_LDFLAGS =  $(AM_LDFLAGS)  @MPI_LIBS@ @BOOST_LDFLAGS@\
                 @BOOST_MPI@ @BOOST_PROGOP@ @BOOST_FS@ @BOOST_REGEX@
cubic_LDADD =   libce.la \
                ../crystal/libcrystal.la \
		../atat/libatat.la \
		../physics/libphysics.la \
		@MPI_ADDLIB@ \
		../print/libprint.la 
tetra_SOURCES = $(cubic_SOURCES)
tetra_CPPFLAGS =  $(AM_CPPFLAGS) -I .. -D_TETRAGONAL_CE_
tetra_LDFLAGS = $(cubic_LDFLAGS)
tetra_LDADD = $(cubic_LDADD) 
regular_SOURCES = rmain.cc
regular_CPPFLAGS =  $(AM_CPPFLAGS) -I .. -D_CUBIC_CE_ @GSL_CFLAGS@
regular_LDFLAGS = $(cubic_LDFLAGS) @GSL_LIBS@ 
regular_LDADD = libce.la \
                ../opt/libopt.la \
                ../crystal/libcrystal.la \
		../atat/libatat.la \
		../physics/libphysics.la \
		../print/libprint.la 
noinst_LTLIBRARIES = libce.la 
libce_la_SOURCES = $(CE_CC)
libce_ladir = $(includedir)/ce
nobase_libce_la_HEADERS = $(CE_H)
libce_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
../crystal/licrystal.la: 
	cd ../crystal; $(MAKE) 

../atat/libatat.la: 
	cd ../atat; $(MAKE) 

../atat/libphysic.la: 
	cd ../physics; $(MAKE) 

../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)

../print/libprint.la:
	cd ../print; $(MAKE)

../opt/libopt.la:
	cd ../opt; $(MAKE)
