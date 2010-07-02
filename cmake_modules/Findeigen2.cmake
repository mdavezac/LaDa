if(eigen2_INCLUDE_DIR)
  set( eigen2_FIND_QUIETLY True)
endif(eigen2_INCLUDE_DIR)
find_path(eigen2_INCLUDE_DIR
  NAMES
  Core
  LU
  Geometry
  Cholesky
  PATHS
  $ENV{eigen2_INCLUDE_DIR}
  $ENV{HOME}/usr/include/
  $ENV{HOME}/hopper/include/
  /usr/include
  /usr/local/include
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  eigen2
)


IF (eigen2_INCLUDE_DIR)
   SET(EIGEN2_FOUND TRUE)
   SET(eigen2_INCLUDE_DIR ${eigen2_INCLUDE_DIR}/eigen2)
ENDIF (eigen2_INCLUDE_DIR)


IF (EIGEN2_FOUND)
   IF (NOT eigen2_FIND_QUIETLY)
      MESSAGE(STATUS "Found eigen2: ${eigen2_INCLUDE_DIR}")
   ENDIF (NOT eigen2_FIND_QUIETLY)
ELSE (EIGEN2_FOUND)
   IF (EIGEN2_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find eigen2")
   ENDIF (EIGEN2_FIND_REQUIRED)
ENDIF (EIGEN2_FOUND)

