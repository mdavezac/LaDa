include ../Files
AM_CPPFLAGS = @COMPILER_CXXFLAGS@  @CPP_DEBUG_FLAGS@
AM_LDFLAGS = 
noinst_LTLIBRARIES = libatat.la
lib_LTLIBRARIES = libatat.py.la
libatat_la_SOURCES = $(ATAT_CC)
libatat_ladir = $(includedir)/atat
nobase_libatat_la_HEADERS = $(ATAT_H)
libatat_py_la_CPPFLAGS = $(AM_CPPFLAGS)  $(AM_CXXFLAGS) @BOOST_CXXFLAGS@ \
       	                 -I@PYTHON_INCLUDE_DIR@
libatat_py_la_SOURCES = $(ATAT_PYTHON_CC)
libatat_py_la_LDFLAGS = $(AM_LDFLAGS)
libatat_py_la_LIBADD = @BOOST_LDFLAGS@ @BOOST_PYTHON@
# libatat_py_la_HEADERS = $(ATAT_PY_HPP)
