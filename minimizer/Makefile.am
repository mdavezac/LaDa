#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @CPP_DEBUG_FLAGS@ @MPI_CFLAGS@ \
	      @BOOST_CXXFLAGS@  @EIGEN_CXXFLAGS@ @TINYXML_CFLAGS@ @MINUIT_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@
noinst_LTLIBRARIES = libminimizer.la
libminimizer_la_SOURCES =  $(minimizer_f90) $(minimizer_cc)
libminimizer_la_CPPFLAGS = $(AM_CPPFLAGS)
libminimizer_ladir = $(includedir)/LaDa/minimizer
libminimizer_la_HEADERS = $(minimizer_h)
SUBDIRS = python
