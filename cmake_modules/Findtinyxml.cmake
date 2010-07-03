if(tinyxml_INCLUDE_DIRS AND tinyxml_LIBRARY)
  set( tinyxml_FIND_QUIETLY True)
endif(tinyxml_INCLUDE_DIRS AND tinyxml_LIBRARY)

find_path(tinyxml_INCLUDE_DIRS
  NAMES
  tinystr.h
  tinyxml.h
  PATHS
  $ENV{tinyxml_INCLUDE_DIRS}
  $ENV{HOME}/usr/include/
  $ENV{HOME}/hopper/include/
  /usr/include
  /usr/local/include
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  tinyxml
)

FIND_LIBRARY(tinyxml_LIBRARY
  tinyxml
  PATH
  $ENV{tinyxml_LIBRARY_DIR}
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
FIND_PACKAGE_HANDLE_STANDARD_ARGS(tinyxml DEFAULT_MSG tinyxml_LIBRARY tinyxml_INCLUDE_DIRS)

IF(tinyxml_INCLUDE_DIRS AND tinyxml_LIBRARY)
  SET(tinyxml_FOUND TRUE)
ENDIF(tinyxml_INCLUDE_DIRS AND tinyxml_LIBRARY)
