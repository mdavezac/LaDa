#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS = -I $(top_srcdir) @CPP_DEBUG_FLAGS@ @MPI_CFLAGS@ \
	      @BOOST_CXXFLAGS@ @TINYXML_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@ @FCLIBS@
noinst_LTLIBRARIES = libopt.la libpythonconvenience.la 
libpythonconvenience_ladir = $(includedir)/LaDa/opt
libpythonconvenience_la_SOURCES = $(OPT_PYTHON_CONVENIENCE_CC)
libpythonconvenience_la_HEADERS = $(OPT_PYTHON_CONVENIENCE_H)
libpythonconvenience_la_CPPFLAGS = $(AM_CPPFLAGS)
libopt_la_SOURCES = $(OPT_F90) $(OPT_CC) 
libopt_la_CPPFLAGS = $(AM_CPPFLAGS)
libopt_la_LDFLAGS = $(AM_LDFLAGS) @FCLIBS@
libopt_la_LIBADD = libpythonconvenience.la
libopt_ladir = $(includedir)/LaDa/opt
nobase_libopt_la_HEADERS = $(OPT_H)
SUBDIRS = python
