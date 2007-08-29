AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @CPP_DEBUG_FLAGS@ @MPICH_CFLAGS@
noinst_LTLIBRARIES = libprint.la
libprint_ladir = print
libprint_la_SOURCES = manip.cc xmg.cc stdout.cc
nobase_libprint_la_HEADERS = manip.h xmg.h stdout.h
libprint_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
