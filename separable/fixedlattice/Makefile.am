#
# Version: $Id$
#
include $(top_srcdir)/Files

bin_PROGRAMS = mixed mixedtetra # tetra cubic 

AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ \
	      @BOOST_CXXFLAGS@ -D __DOHALFHALF__
AM_LDFLAGS = @TINYXML_LIBS@ @BOOST_LDFLAGS@ @BOOST_PROGOP@ \
             @BOOST_FS@ @BOOST_REGEX@ @BOOST_SYSTEM@ 

# cubic_SOURCES = main.cc  
# cubic_CPPFLAGS =  $(AM_CPPFLAGS)  -D _CUBIC_CE_
# cubic_LDFLAGS =  $(AM_LDFLAGS)
# cubic_LDADD = libfxseparable.la \
# 	      ../libseparable.la \
# 	      ../../crystal/libcrystal.la \
# 	      ../../atat/libatat.la \
# 	      ../../physics/libphysics.la \
# 	      ../../print/libprint.la \
# 	      ../../opt/libopt.la  
# 
# tetra_SOURCES = $(cubic_SOURCES)
# tetra_CPPFLAGS =  $(AM_CPPFLAGS) -D_TETRAGONAL_CE_
# tetra_LDFLAGS = $(cubic_LDFLAGS)
# tetra_LDADD = $(cubic_LDADD) 

mixed_SOURCES = mxmain.cc
mixed_CPPFLAGS =  $(AM_CPPFLAGS) -D_CUBIC_CE_
mixed_LDFLAGS = $(AM_LDFLAGS)
mixed_LDADD = libfxseparable.la \
              ../libseparable.la \
              ../../ce/libce.la \
              ../../crystal/libcrystal.la \
              ../../atat/libatat.la \
              ../../physics/libphysics.la \
              ../../print/libprint.la \
              ../../opt/libopt.la  

mixedtetra_SOURCES = mxmain.cc
mixedtetra_CPPFLAGS =  $(AM_CPPFLAGS) -D_TETRAGONAL_CE_
mixedtetra_LDFLAGS = $(mixed_LDFLAGS)
mixedtetra_LDADD = $(mixed_LDADD) 

noinst_LTLIBRARIES = libfxseparable.la 
libfxseparable_la_SOURCES = $(SEPARABLES_FXL_CC)
libfxseparable_ladir = $(includedir)/LaDa/separable/fixedlattice
nobase_libfxseparable_la_HEADERS = $(SEPARABLES_FXL_H) $(SEPARABLES_FXL_IMPL_H)
libfxseparable_la_CPPFLAGS = $(AM_CPPFLAGS) 

../libseparable.la:
	cd ../; $(MAKE) libseparable.la
../../atat/libatat.la: 
	cd ../../atat; $(MAKE) libatat.la
../../crystal/libcrystal.la: 
	cd ../../crystal; $(MAKE) libcrystal.la
../../physics/libphysic.la: 
	cd ../../physics; $(MAKE) libphysics.la
../../opt/libopt.la: 
	cd ../../opt; $(MAKE) libopt.la
../../print/libprint.la: 
	cd ../../print; $(MAKE) libprint.la
../../ce/libce.la: 
	cd ../../ce; $(MAKE) libce.la
../../mpi/libmpiobject.la:
	cd ../../mpi; $(MAKE) libmpiobject.la
