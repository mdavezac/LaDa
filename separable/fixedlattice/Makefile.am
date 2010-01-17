#
# Version: $Id$
#
include $(top_srcdir)/Files

bin_PROGRAMS = mixed mixedtetra # tetra cubic 

AM_CPPFLAGS = -I $(top_srcdir) @TINYXML_CFLAGS@ @CPP_DEBUG_FLAGS@ \
	      @BOOST_CXXFLAGS@ -D __DOHALFHALF__
AM_LDFLAGS = @TINYXML_LIBS@ @BOOST_LDFLAGS@ @BOOST_PROGOP@ \
             @BOOST_FS@ @BOOST_REGEX@ @BOOST_SYSTEM@ @BOOST_MPI@

# cubic_SOURCES = main.cc  
# cubic_CPPFLAGS =  $(AM_CPPFLAGS)  -D _CUBIC_CE_
# cubic_LDFLAGS =  $(AM_LDFLAGS)
# cubic_LDADD = libfxseparable.la \
# 	      ../libseparable.la \
# 	      ../../crystal/libcrystal.la \
# 	      ../../math/libmath.la \
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
              $(top_builddir)/separable/libseparable.la \
              $(top_builddir)/ce/libce.la \
              $(top_builddir)/crystal/libcrystal.la \
              $(top_builddir)/math/libmath.la \
              $(top_builddir)/physics/libphysics.la \
              $(top_builddir)/print/libprint.la \
              $(top_builddir)/opt/libopt.la  \
              $(top_builddir)/mpi/libmpiobject.la  

mixedtetra_SOURCES = mxmain.cc
mixedtetra_CPPFLAGS =  $(AM_CPPFLAGS) -D_TETRAGONAL_CE_
mixedtetra_LDFLAGS = $(mixed_LDFLAGS)
mixedtetra_LDADD = $(mixed_LDADD) 

noinst_LTLIBRARIES = libfxseparable.la 
libfxseparable_la_SOURCES = $(separable_fxl_cc)
libfxseparable_ladir = $(includedir)/LaDa/separable/fixedlattice
nobase_libfxseparable_la_HEADERS = $(separable_fxl_h) $(separable_fxl_impl_h)
libfxseparable_la_CPPFLAGS = $(AM_CPPFLAGS) 
