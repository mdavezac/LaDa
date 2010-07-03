if(eigen2_INCLUDE_DIR)
  set( eigen2_FIND_QUIETLY True)
endif(eigen2_INCLUDE_DIR)

find_path(eigen2_INCLUDE_DIR
  NAMES
  Eigen/Core
  Eigen/LU
  Eigen/Geometry
  Eigen/Cholesky
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


INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(eigen2 DEFAULT_MSG eigen2_INCLUDE_DIR)
