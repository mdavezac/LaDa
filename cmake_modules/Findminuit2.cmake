if(minuit2_INCLUDE_DIRS AND minuit2_LIBRARY)
  set( minuit2_FIND_QUIETLY True)
endif(minuit2_INCLUDE_DIRS AND minuit2_LIBRARY)

find_path(minuit2_INCLUDE_DIRS
  NAMES
  MatrixInverse.h
  Minuit2Minimizer.h
  PATHS
  $ENV{minuit2_INCLUDE_DIRS}
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  Minuit2
)

FIND_LIBRARY(minuit2_LIBRARY
  Minuit2
  PATH
  $ENV{minuit2_LIBRARY_DIR}
)


INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(minuit2 DEFAULT_MSG minuit2_LIBRARY minuit2_INCLUDE_DIRS)
