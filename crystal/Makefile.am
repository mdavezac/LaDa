#
# Version: $Id$
#
include $(top_srcdir)/Files
# if DOPYTHON
#   pythonlib = libcrystalpy.la
# endif
AM_CPPFLAGS = @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
	      @BOOST_CXXFLAGS@ -I $(top_srcdir) 
AM_LDFLAGS = @BOOST_REGEX@ @BOOST_FS@ @BOOST_SYSTEM@ \
             @BOOST_MPI@ @BOOST_SERIALIZATION@

noinst_LTLIBRARIES = libcrystal.la # $(pythonlib)
libcrystal_la_SOURCES = $(crystal_cc)
libcrystal_ladir = $(includedir)/LaDa/crystal
nobase_libcrystal_la_HEADERS = $(crystal_h) 
libcrystal_la_CPPFLAGS = $(AM_CPPFLAGS) 
SUBDIRS = python

# if DOALLOYLAYERSOPT
# bin_PROGRAMS = alloylayers
# alloylayers_SOURCES = main.cc
# alloylayers_CPPFLAGS =  $(AM_CPPFLAGS)
# alloylayers_LDFLAGS =  $(AM_LDFLAGS) @BOOST_PROGOP@ 
# alloylayers_LDADD = libcrystal.la \
#	      $(top_builddir)/physics/libphysics.la \
#	      $(top_builddir)/math/libmath.la \
#	      $(top_builddir)/opt/libopt.la \
#	      $(top_builddir)/mpi/libmpiobject.la \
#	      $(top_builddir)/print/libprint.la \
#	      @BOOST_MPI@ @BOOST_SERIALIZATION@ @BOOST_REGEX@
