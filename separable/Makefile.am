include ../Files
bin_PROGRAMS = separables
AM_CPPFLAGS = @COMPILER_CXXFLAGS@ @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ \
	      @BOOST_CXXFLAGS@ @GSL_CFLAGS@
AM_LDFLAGS = @TINYXML_LIBS@ @BOOST_LDFLAGS@ @BOOST_PROGOP@ \
             @BOOST_FS@ @GSL_LIBS@
separables_SOURCES = main.cc 
separables_CPPFLAGS =  $(AM_CPPFLAGS) -I .. 
separables_LDFLAGS =  $(AM_LDFLAGS)
separables_LDADD = libseparable.la \
                ../lamarck/liblamarck.la \
		../atat/libatat.la \
		../opt/libopt.la \
		../physics/libphysics.la 
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
