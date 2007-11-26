AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @CPP_DEBUG_FLAGS@ @MPICH_CFLAGS@
noinst_LTLIBRARIES = libprint.la
libprint_ladir = $(includedir)/print
libprint_la_SOURCES = manip.cc xmg.cc stdout.cc base.cc 
nobase_libprint_la_HEADERS = manip.h base.h xmg.h stdout.h operations.h
libprint_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 
