include ../Files

bin_PROGRAMS = tetra cubic

AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ \
	      @BOOST_CXXFLAGS@ @GSL_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@ @BOOST_LDFLAGS@ @BOOST_PROGOP@ \
             @BOOST_FS@ @GSL_LIBS@ @BOOST_REGEX@

cubic_SOURCES = main.cc 
cubic_CPPFLAGS =  $(AM_CPPFLAGS) -I ..  -D _CUBIC_CE_
cubic_LDFLAGS =  $(AM_LDFLAGS)
cubic_LDADD = libseparable.la \
              ../lamarck/liblamarck.la \
              ../atat/libatat.la \
              ../opt/libopt.la \
              ../physics/libphysics.la \
              ../print/libprint.la 

tetra_SOURCES = $(cubic_SOURCES)
tetra_CPPFLAGS =  $(AM_CPPFLAGS) -I .. -D_TETRAGONAL_CE_
tetra_LDFLAGS = $(cubic_LDFLAGS)
tetra_LDADD = $(cubic_LDADD) 

noinst_LTLIBRARIES = libseparable.la 
libseparable_la_SOURCES = $(SEPARABLE_CC)
libseparable_ladir = $(includedir)/separable
nobase_libseparable_la_HEADERS = $(SEPARABLE_H) $(SEPARABLE_IMPL_H)
libseparable_la_CPPFLAGS = -I ..  $(AM_CPPFLAGS) 

../atat/libatat.la: 
	cd ../atat; $(MAKE) 
../lamarck/liblamarck.la: 
	cd ../lamarck; $(MAKE) 
../physics/libphysic.la: 
	cd ../physics; $(MAKE) 
../opt/libopt.la: 
	cd ../opt; $(MAKE) 
../opt/libprint.la: 
	cd ../print; $(MAKE) 
