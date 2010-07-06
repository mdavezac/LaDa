# - Check for the presence of PYUBLAS
#
# The following variables are set when NUMPY is found:
#  PYUBLAS_FOUND              = Set to true, if all components of NUMPY have been
#                               found.
#  PYUBLAS_INCLUDES           = Include path for the header files of NUMPY

# standard required packages
find_package(PythonLibs REQUIRED)
find_package(PythonInterp REQUIRED)
find_package(numpy REQUIRED)

## -----------------------------------------------------------------------------
## Let's assume the Python executable is smarter about finding NumPy than we
## are, and try asking it before searching ourselves.
## This is necessary to e.g. pick up the MacPorts NumPy installation, which
## ends up in /opt/local/Library/Frameworks/Python.framework ...

execute_process (
  COMMAND ${PYTHON_EXECUTABLE} -c "import pyublas, os; print os.path.dirname(pyublas.__file__)"
  OUTPUT_VARIABLE pyublas_path
  )
if (pyublas_path)
  string (STRIP ${pyublas_path} pyublas_search_path)
else (pyublas_path)
  set (pyublas_search_path ${lib_locations})
endif (pyublas_path)

## -----------------------------------------------------------------------------
## Check for the header files

find_path (PYUBLAS_INCLUDE_DIRS pyublas/numpy.hpp
  PATHS
  ${pyublas_search_path}
  PATH_SUFFIXES
  ../include
  )


## -----------------------------------------------------------------------------
## Done

unset(pyublas_path CACHE)
set(PYUBLAS_INCLUDE_DIRS ${PYUBLAS_INCLUDE_DIRS} CACHE PATH "Include directory of PyUblas.")

