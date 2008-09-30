#
# Version: $Id$
#
include ../Files
AM_CPPFLAGS = @CPP_DEBUG_FLAGS@ @MPI_CFLAGS@ @BOOST_CXXFLAGS@
AM_LDFLAGS = @BOOST_LDFLAGS@
noinst_LTLIBRARIES = libmpiobject.la
libmpiobject_la_SOURCES = $(MPI_CC)
libmpiobject_la_HEADERS = $(MPI_H) 
libmpiobject_ladir = $(includedir)/LaDa/mpi
libmpiobject_la_LDFLAGS = $(AM_LDFLAGS) @MPI_LIBS@
libmpiobject_la_CPPFLAGS = $(AM_CPPFLAGS) $(MPI_CFLAGS)



bin_PROGRAMS = mpitest
mpitest_SOURCES = main.cc 
mpitest_LDFLAGS =  $(AM_LDFLAGS) @MPI_LIBS@ @BOOST_MPI@
mpitest_LDADD = @MPI_ADDLIB@ ../print/libprint.la
../print/libprint.la: 
	cd ../print; $(MAKE)
