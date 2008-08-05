#
# Version: $Id$
#
include ../Files
if DOPYTHON
  pythonlib = libatatpy.la
endif
AM_CPPFLAGS = @COMPILER_CXXFLAGS@  @CPP_DEBUG_FLAGS@ @BOOST_CXXFLAGS@ 
AM_LDFLAGS = 
noinst_LTLIBRARIES = libatat.la $(pythonlib)
libatat_la_SOURCES = $(ATAT_CC)
libatat_ladir = $(includedir)/LaDa/atat
nobase_libatat_la_HEADERS = $(ATAT_H)
if DOPYTHON
libatatpy_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) \
                        -I@PYTHON_INCLUDE_DIR@
libatatpy_la_SOURCES = $(ATAT_PYTHON_CC)
libatatpy_la_LDFLAGS = $(AM_LDFLAGS) -Bshared 
libatatpy_la_LIBADD = @BOOST_LDFLAGS@ @BOOST_PYTHON@
endif
