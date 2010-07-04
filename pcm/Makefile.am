#
# Version: $Id$
#
include $(top_srcdir)/Files
# if DOPYTHON
#   pythonlib = libcrystalpy.la
# endif
AM_CPPFLAGS = @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@\
      	@BOOST_CXXFLAGS@ -I $(top_srcdir) @GSL_CFLAGS@ @MINUIT_CFLAGS@
AM_LDFLAGS = @BOOST_FS@ @BOOST_SYSTEM@ @GSL_LIBS@ @TINYXML_LIBS@ @FCLIBS@ @MINUIT_LIBS@

lib_LTLIBRARIES = librelaxer.la
noinst_LTLIBRARIES = libmodels.la

libmodels_la_SOURCES = $(models_clj_cc) $(models_clj_f90)
libmodels_la_CPPFLAGS = $(AM_CPPFLAGS) 
libmodels_ladir = $(includedir)/LaDa/models
nobase_libmodels_la_HEADERS = $(models_clj_h) 


EXTRA_DIST = $(fc_models_clj_h) 
librelaxer_la_SOURCES = $(fc_models_clj_cc)
librelaxer_la_CPPFLAGS = $(AM_CPPFLAGS) 
librelaxer_la_LIBADD = $(top_builddir)/crystal/libcrystal.la \
                       $(top_builddir)/math/libmath.la \
                       $(top_builddir)/minimizer/libminimizer.la \
                       $(top_builddir)/opt/libopt.la \
                       $(top_builddir)/physics/libphysics.la

# bin_PROGRAMS = clj coollj
# clj_SOURCES =  CLJ.f90 force_lj_ewald.f90
# clj_LDADD = libmodels.la librelaxer.la  @BOOST_PROGOP@ 
# coollj_SOURCES =  main.cc
# coollj_LDADD = @BOOST_PROGOP@ libmodels.la \
#	 $(top_builddir)/crystal/libcrystal.la \
#	 $(top_builddir)/math/libmath.la \
#	 $(top_builddir)/minimizer/libminimizer.la \
#	 $(top_builddir)/opt/libopt.la \
#	 $(top_builddir)/physics/libphysics.la


SUBDIRS = python
