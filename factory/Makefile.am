#
# Version: $Id: Makefile.am 851 2008-11-13 01:57:17Z davezac $
#
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @CPP_DEBUG_FLAGS@ @MPI_CFLAGS@ \
	      @BOOST_CXXFLAGS@ @TINYXML_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@ @BOOST_LDFLAGS@ @BOOST_SIGNAL@ 
noinst_LTLIBRARIES = libfactory.la
libfactory_ladir = $(includedir)/LaDa/factory
libfactory_la_HEADERS = $(FACTORY_H)
libfactory_la_SOURCES = 
