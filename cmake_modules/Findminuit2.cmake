if(minuit2_INCLUDE_DIRS)
  set( minuit2_FIND_QUIETLY True)
endif(minuit2_INCLUDE_DIRS)

find_path(minuit2_INCLUDE_DIRS
  NAMES
  MatrixInverse.h
  Minuit2Minimizer.h
  PATHS
  $ENV{minuit2_INCLUDE_DIRS}
  $ENV{HOME}/usr/include/
  $ENV{HOME}/hopper/include/
  /usr/include
  /usr/local/include
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  Minuit2
)

FIND_LIBRARY(minuit2_LIBRARY
  Minuit2
  PATH
  $ENV{minuit2_LIBRARY_DIR}
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
FIND_PACKAGE_HANDLE_STANDARD_ARGS(minuit2 DEFAULT_MSG minuit2_LIBRARY minuit2_INCLUDE_DIRS)

IF(minuit2_INCLUDE_DIRS AND minuit2_LIBRARY)
  SET(minuit2_FOUND TRUE)
ENDIF(minuit2_INCLUDE_DIRS AND minuit2_LIBRARY)
