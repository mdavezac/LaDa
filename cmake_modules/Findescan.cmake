if(ESCAN_INCLUDE_DIRS AND ESCAN_LIBRARY AND GENPOT_LIBRARY)
  set( ESCAN_FIND_QUIETLY TRUE)
endif(ESCAN_INCLUDE_DIRS AND ESCAN_LIBRARY AND GENPOT_LIBRARY)

find_path(_ESCAN_INCLUDE_DIRS
  NAMES
  coulomb_4pair_api.mod
  PATHS
  $ENV{ESCAN_INCLUDE_DIRS}
  $ENV{HOME}/usr/include/
  $ENV{HOME}/hopper/include/
  /usr/include
  /usr/local/include
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  escan nanopse
)
if(NOT _ESCAN_INCLUDE_DIRS)
  find_path(_ESCAN_INCLUDE_DIRS
    NAMES
    COULOMB_4PAIR_API.mod
    PATHS
    $ENV{ESCAN_INCLUDE_DIRS}
    $ENV{HOME}/usr/include/
    $ENV{HOME}/hopper/include/
    /usr/include
    /usr/local/include
    ${INCLUDE_INSTALL_DIR}
    PATH_SUFFIXES
    escan nanopse
  )
endif(NOT _ESCAN_INCLUDE_DIRS)


FIND_LIBRARY(_ESCAN_LIBRARY
  pescan
  PATH
  $ENV{ESCAN_LIBRARY_DIR}
  $ENV{HOME}/usr/lib64/
  $ENV{HOME}/usr/lib/
  $ENV{HOME}/hopper/lib64
  $ENV{HOME}/hopper/lib
  /usr/lib64
  /usr/lib
  /usr/local/lib64
  /usr/local/lib
)
FIND_LIBRARY(_GENPOT_LIBRARY
  genpot
  PATH
  $ENV{ESCAN_LIBRARY_DIR}
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
FIND_PACKAGE_HANDLE_STANDARD_ARGS(escan DEFAULT_MSG _ESCAN_LIBRARY _GENPOT_LIBRARY _ESCAN_INCLUDE_DIRS)
IF(_ESCAN_INCLUDE_DIRS AND _ESCAN_LIBRARY AND _GENPOT_LIBRARY)
  set(ESCAN_INCLUDE_DIRS ${_ESCAN_INCLUDE_DIRS} CACHE PATH "Path to nanopse include directory.")
  set(ESCAN_LIBRARY ${_ESCAN_LIBRARY} CACHE PATH "Path to escan library.")
  set(GENPOT_LIBRARY ${_GENPOT_LIBRARY} CACHE PATH "Path to genpot library.")
  SET(ESCAN_FOUND TRUE)
  unset(_ESCAN_INCLUDE_DIRS CACHE)
  unset(_GENPOT_LIBRARY CACHE)
  unset(_ESCAN_LIBRARY CACHE)
ENDIF(_ESCAN_INCLUDE_DIRS AND _ESCAN_LIBRARY AND _GENPOT_LIBRARY)
