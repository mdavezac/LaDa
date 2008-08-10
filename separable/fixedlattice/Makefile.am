#
# Version: $Id$
#
include ../../Files

bin_PROGRAMS = tetra cubic

AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ \
	      @BOOST_CXXFLAGS@ -D __DOHALFHALF__
AM_LDFLAGS = @TINYXML_LIBS@ @BOOST_LDFLAGS@ @BOOST_PROGOP@ \
             @BOOST_FS@ @BOOST_REGEX@ @BOOST_SYSTEM@

cubic_SOURCES = main.cc  
cubic_CPPFLAGS =  $(AM_CPPFLAGS) -I ..  -D _CUBIC_CE_
cubic_LDFLAGS =  $(AM_LDFLAGS)
cubic_LDADD = libseparable.la \
              ../../crystal/libcrystal.la \
              ../../atat/libatat.la \
              ../../opt/libopt.la \
              ../../physics/libphysics.la \
              ../../print/libprint.la 

tetra_SOURCES = $(cubic_SOURCES)
tetra_CPPFLAGS =  $(AM_CPPFLAGS) -I .. -D_TETRAGONAL_CE_
tetra_LDFLAGS = $(cubic_LDFLAGS)
tetra_LDADD = $(cubic_LDADD) 

noinst_LTLIBRARIES = libseparable.la 
libseparable_la_SOURCES = $(SEPARABLES_FXL_CC)
libseparable_ladir = $(includedir)/LaDa/separable/fixedlattice
nobase_libseparable_la_HEADERS = $(SEPARABLES_FXL_H) $(SEPARABLES_FXL_IMPL_H)
libseparable_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 

SUBDIR = fixedlattice

../../atat/libatat.la: 
	cd ../atat; $(MAKE) 
../../crystal/libcrystal.la: 
	cd ../crystal; $(MAKE) 
../../physics/libphysic.la: 
	cd ../physics; $(MAKE) 
../../opt/libopt.la: 
	cd ../opt; $(MAKE) 
../../opt/libprint.la: 
	cd ../print; $(MAKE) 
