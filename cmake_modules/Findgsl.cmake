if(gsl_INCLUDE_DIRS)
  set( gsl_FIND_QUIETLY True)
endif(gsl_INCLUDE_DIRS)

find_path(gsl_INCLUDE_DIRS
  NAMES
  gsl_vector.h
  gsl_matrix.h
  PATHS
  $ENV{gsl_INCLUDE_DIRS}
  $ENV{HOME}/usr/include/
  $ENV{HOME}/hopper/include/
  /usr/include
  /usr/local/include
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  gsl
)

FIND_LIBRARY(gsl_LIBRARY
  gsl
  gslcblas
  PATH
  $ENV{gsl_LIBRARY_DIR}
  $ENV{HOME}/usr/lib64/
  $ENV{HOME}/usr/lib/
  $ENV{HOME}/hopper/lib64
  $ENV{HOME}/hopper/lib
  /usr/lib64
  /usr/lib
  /usr/local/lib64
  /usr/local/lib
)


INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(gsl DEFAULT_MSG gsl_LIBRARY gsl_INCLUDE_DIRS)

IF(gsl_INCLUDE_DIRS AND gsl_LIBRARY)
  SET(gsl_FOUND TRUE)
ENDIF(gsl_INCLUDE_DIRS AND gsl_LIBRARY)
