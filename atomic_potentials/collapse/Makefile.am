#
# Version: $Id$
#
include $(top_srcdir)/Files

AM_CPPFLAGS = @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ @BOOST_CXXFLAGS@ -I $(top_srcdir) 
AM_LDFLAGS = 

noinst_LTLIBRARIES = libcollapse.la
libcollapse_la_SOURCES = $(ap_collapse_cc)
libcollapse_ladir = $(includedir)/LaDa/atomic_potentials/collapse
nobase_libcollapse_la_HEADERS = $(ap_collapse_h)
libcollapse_la_CPPFLAGS = $(AM_CPPFLAGS) 
