if(escan_INCLUDE_DIRS AND escan_LIBRARY)
  set( escan_FIND_QUIETLY True)
endif(escan_INCLUDE_DIRS AND escan_LIBRARY)

find_path(escan_INCLUDE_DIRS
  NAMES
  coulomb_4pair_api.mod
  PATHS
  $ENV{escan_INCLUDE_DIRS}
  $ENV{HOME}/usr/include/
  $ENV{HOME}/hopper/include/
  /usr/include
  /usr/local/include
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  escan nanopse
)
if(NOT escan_INCLUDE_DIRS)
  find_path(escan_INCLUDE_DIRS
    NAMES
    COULOMB_4PAIR_API.mod
    PATHS
    $ENV{escan_INCLUDE_DIRS}
    $ENV{HOME}/usr/include/
    $ENV{HOME}/hopper/include/
    /usr/include
    /usr/local/include
    ${INCLUDE_INSTALL_DIR}
    PATH_SUFFIXES
    escan nanopse
  )
endif(NOT escan_INCLUDE_DIRS)


FIND_LIBRARY(escan_LIBRARY
  pescan
  genpot
  PATH
  $ENV{escan_LIBRARY_DIR}
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
FIND_PACKAGE_HANDLE_STANDARD_ARGS(escan DEFAULT_MSG escan_LIBRARY escan_INCLUDE_DIRS)

IF(escan_INCLUDE_DIRS AND escan_LIBRARY)
  SET(escan_FOUND TRUE)
ENDIF(escan_INCLUDE_DIRS AND escan_LIBRARY)
