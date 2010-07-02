#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS = @CPP_DEBUG_FLAGS@ @MPI_CFLAGS@ @BOOST_CXXFLAGS@ -I$(top_srcdir)
AM_LDFLAGS = @BOOST_LDFLAGS@
noinst_LTLIBRARIES = libmpiobject.la
libmpiobject_la_SOURCES = $(mpi_cc)
libmpiobject_la_HEADERS = $(mpi_h) 
libmpiobject_ladir = $(includedir)/LaDa/mpi
libmpiobject_la_LDFLAGS = $(AM_LDFLAGS) 
libmpiobject_la_CPPFLAGS = $(AM_CPPFLAGS) $(MPI_CFLAGS)

