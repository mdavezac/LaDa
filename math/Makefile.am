#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS = -I$(top_srcdir)  @CPP_DEBUG_FLAGS@ @BOOST_CXXFLAGS@ 
AM_LDFLAGS = 
noinst_LTLIBRARIES = libmath.la 
libmath_la_SOURCES = $(math_cc) $(math_f90)
libmath_ladir = $(includedir)/LaDa/math
nobase_libmath_la_HEADERS = $(math_h)
# SUBDIRS = python
