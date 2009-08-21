#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS = -I$(top_srcdir)  @CPP_DEBUG_FLAGS@ @BOOST_CXXFLAGS@ 
AM_LDFLAGS = 
noinst_LTLIBRARIES = libatat.la 
libatat_la_SOURCES = $(atat_cc)
libatat_ladir = $(includedir)/LaDa/atat
nobase_libatat_la_HEADERS = $(atat_h)
SUBDIRS = python
