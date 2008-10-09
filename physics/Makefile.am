#
# Version: $Id$
#
include $(top_srcdir)/Files
AM_CPPFLAGS =  @CPP_DEBUG_FLAGS@  -I $(top_srcdir)
noinst_LTLIBRARIES = libphysics.la
libphysics_la_SOURCES = $(PHYSICS_CC)
libphysics_la_CPPFLAGS = $(AM_CPPFLAGS) $(AM_CXXFLAGS) @BOOST_CXXFLAGS@
libphysics_ladir = $(includedir)/LaDa/physics
nobase_libphysics_la_HEADERS = $(PHYSICS_H)
