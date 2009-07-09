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

libmodels_la_SOURCES = $(MODELS_CLJ_CC) $(MODELS_CLJ_F90)
libmodels_la_CPPFLAGS = $(AM_CPPFLAGS) 
libmodels_ladir = $(includedir)/LaDa/models
nobase_libmodels_la_HEADERS = $(MODELS_CLJ_H) 


librelaxer_la_SOURCES = $(FC_MODELS_CLJ_CC)
librelaxer_la_CPPFLAGS = $(AM_CPPFLAGS) 
librelaxer_la_LIBADD = ../crystal/libcrystal.la \
                       ../atat/libatat.la \
                       ../minimizer/libminimizer.la \
                       ../opt/libopt.la \
                       ../print/libprint.la \
                       ../physics/libphysics.la

bin_PROGRAMS = # clj coollj
clj_SOURCES =  CLJ.f90 force_lj_ewald.f90
clj_LDADD = libmodels.la librelaxer.la  @BOOST_PROGOP@ 
coollj_SOURCES =  main.cc
coollj_LDADD = @BOOST_PROGOP@ libmodels.la \
               ../crystal/libcrystal.la \
               ../atat/libatat.la \
               ../minimizer/libminimizer.la \
               ../opt/libopt.la \
               ../print/libprint.la \
               ../physics/libphysics.la


../minimizer/libminimizer.la: 
	cd ../opt; $(MAKE) libminimizer.la

../opt/libopt.la: 
	cd ../opt; $(MAKE) libopt.la

../physics/libphysics.la: 
	cd ../physics; $(MAKE) libphysics.la

../print/libprint.la: 
	cd ../print; $(MAKE) libprint.la

../atat/libatat.la: 
	cd ../atat; $(MAKE) libatat.la

../crystal/libcrystal.la:
	cd ../crystal; $(MAKE) libcrystal.la

SUBDIRS = python
