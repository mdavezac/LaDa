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
librelaxer_la_SOURCES = $(MODELS_CLJ_CC)
librelaxer_ladir = $(includedir)/LaDa/lennard-jhones
librelaxer_la_HEADERS = $(MODELS_CLJ_H) 
librelaxer_la_CPPFLAGS = $(AM_CPPFLAGS) 
librelaxer_la_LIBADD = ../crystal/libcrystal.la \
      		 ../atat/libatat.la \
      		 ../minimizer/libminimizer.la \
      		 ../opt/libopt.la \
      		 ../print/libprint.la \
      		 ../physics/libphysics.la

bin_PROGRAMS = clj
clj_SOURCES =  ep_com.f90 CLJ.f90 ewaldf.f90  force_lj_ewald.f90
clj_LDADD = librelaxer.la  @BOOST_PROGOP@ 


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
