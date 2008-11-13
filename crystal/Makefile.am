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
libcrystal_la_SOURCES = $(CRYSTAL_CC)
libcrystal_ladir = $(includedir)/LaDa/crystal
nobase_libcrystal_la_HEADERS = $(CRYSTAL_H) 
libcrystal_la_CPPFLAGS = $(AM_CPPFLAGS) 

if DOALLOYLAYERSOPT
bin_PROGRAMS = alloylayers
alloylayers_SOURCES = main.cc
alloylayers_CPPFLAGS =  $(AM_CPPFLAGS)
alloylayers_LDFLAGS =  $(AM_LDFLAGS) @BOOST_PROGOP@ 
alloylayers_LDADD = libcrystal.la \
		    ../physics/libphysics.la \
		    ../atat/libatat.la \
		    ../opt/libopt.la \
                    ../mpi/libmpiobject.la
endif

../atat/libatat.la: 
	cd ../atat; $(MAKE) 

../atat/libphysic.la: 
	cd ../physics; $(MAKE) 

../mpi/libmpiobject.la:
	cd ../mpi; $(MAKE)

../print/libprint.la:
	cd ../print; $(MAKE)
