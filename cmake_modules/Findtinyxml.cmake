if(tinyxml_INCLUDE_DIRS AND tinyxml_LIBRARY)
  set( tinyxml_FIND_QUIETLY True)
endif(tinyxml_INCLUDE_DIRS AND tinyxml_LIBRARY)

find_path(tinyxml_INCLUDE_DIRS
  NAMES
  tinystr.h
  tinyxml.h
  PATHS
  $ENV{tinyxml_INCLUDE_DIRS}
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  tinyxml
)

FIND_LIBRARY(tinyxml_LIBRARY
  tinyxml
  PATH
  $ENV{tinyxml_LIBRARY_DIR}
)


INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(tinyxml DEFAULT_MSG tinyxml_LIBRARY tinyxml_INCLUDE_DIRS)

IF(tinyxml_INCLUDE_DIRS AND tinyxml_LIBRARY)
  SET(tinyxml_FOUND TRUE)
ENDIF(tinyxml_INCLUDE_DIRS AND tinyxml_LIBRARY)
