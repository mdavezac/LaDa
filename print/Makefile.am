#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @CPP_DEBUG_FLAGS@ @MPI_CFLAGS@ @BOOST_CXXFLAGS@
noinst_LTLIBRARIES = libprint.la
libprint_ladir = $(includedir)/LaDa/print
libprint_la_SOURCES = $(PRINT_CC)
nobase_libprint_la_HEADERS = $(PRINT_H)
libprint_la_CPPFLAGS = -I $(top_srcdir)  $(AM_CPPFLAGS) 
