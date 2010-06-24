#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS = @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
	      @BOOST_CXXFLAGS@ -I $(top_srcdir) @EIGEN_CXXFLAGS@
AM_LDFLAGS = 

noinst_LTLIBRARIES = libenum.la # $(pythonlib)
libenum_la_SOURCES = $(enum_cc)
libenum_ladir = $(includedir)/LaDa/enum
nobase_libenum_la_HEADERS = $(enum_h) 
libenum_la_CPPFLAGS = $(AM_CPPFLAGS) 
SUBDIRS = python

bin_PROGRAMS = enum
enum_SOURCES = main.cc
enum_CPPFLAGS =  $(AM_CPPFLAGS)
enum_LDFLAGS =  $(AM_LDFLAGS) @BOOST_PROGOP@ 
enum_LDADD = libenum.la \
             $(top_builddir)/crystal/libcrystal.la \
             $(top_builddir)/physics/libphysics.la \
             $(top_builddir)/math/libmath.la \
             $(top_builddir)/opt/libopt.la \
             $(top_builddir)/print/libprint.la 
