#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS = @CPP_DEBUG_FLAGS@ @MPI_CFLAGS@ @BOOST_CXXFLAGS@ -I$(top_srcdir)
AM_LDFLAGS = @BOOST_LDFLAGS@
noinst_LTLIBRARIES = libmpiobject.la
libmpiobject_la_SOURCES = $(MPI_CC)
libmpiobject_la_HEADERS = $(MPI_H) 
libmpiobject_ladir = $(includedir)/LaDa/mpi
libmpiobject_la_LDFLAGS = $(AM_LDFLAGS) 
libmpiobject_la_CPPFLAGS = $(AM_CPPFLAGS) $(MPI_CFLAGS)



bin_PROGRAMS = mpitest
mpitest_SOURCES = main.cc 
mpitest_LDFLAGS =  $(AM_LDFLAGS) @BOOST_MPI@ @BOOST_SERIALIZATION@
mpitest_LDADD = @MPI_ADDLIB@ ../print/libprint.la
../print/libprint.la: 
	cd ../print; $(MAKE)
