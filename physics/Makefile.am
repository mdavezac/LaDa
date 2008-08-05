#
# Version: $Id$
#
include ../Files
AM_CPPFLAGS = @COMPILER_CXXFLAGS@  @CPP_DEBUG_FLAGS@ 
noinst_LTLIBRARIES = libphysics.la
libphysics_la_SOURCES = $(PHYSICS_CC)
libphysics_la_CPPFLAGS = $(AM_CPPFLAGS) $(AM_CXXFLAGS)
libphysics_ladir = $(includedir)/LaDa/physics
nobase_libphysics_la_HEADERS = $(PHYSICS_H)
