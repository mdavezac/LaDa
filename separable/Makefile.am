#
# Version: $Id$
#
include $(top_srcdir)/Files

bin_PROGRAMS = tetra cubic

AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ \
	      @BOOST_CXXFLAGS@ -D __DOHALFHALF__
AM_LDFLAGS = @TINYXML_LIBS@ @BOOST_LDFLAGS@ @BOOST_PROGOP@ \
             @BOOST_SYSTEM@ @BOOST_REGEX@ @BOOST_FS@ @BOOST_MPI@

cubic_SOURCES = main.cc  
cubic_CPPFLAGS =  $(AM_CPPFLAGS) -D _CUBIC_CE_
cubic_LDFLAGS =  $(AM_LDFLAGS)
cubic_LDADD = libseparable.la \
              $(top_builddir)/crystal/libcrystal.la \
              $(top_builddir)/math/libmath.la \
              $(top_builddir)/opt/libopt.la \
              $(top_builddir)/mpi/libmpiobject.la \
              $(top_builddir)/physics/libphysics.la 

tetra_SOURCES = $(cubic_SOURCES)
tetra_CPPFLAGS =  $(AM_CPPFLAGS) -D_TETRAGONAL_CE_
tetra_LDFLAGS = $(cubic_LDFLAGS)
tetra_LDADD = $(cubic_LDADD) 

noinst_LTLIBRARIES = libseparable.la 
libseparable_la_SOURCES = $(separable_cc)
libseparable_ladir = $(includedir)/LaDa/separable
nobase_libseparable_la_HEADERS = $(separable_h) $(separable_impl_h)
libseparable_la_CPPFLAGS = $(AM_CPPFLAGS) 

SUBDIRS = fixedlattice
