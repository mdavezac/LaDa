if(VASP_LIBRARY)
  set(VASP_FIND_QUIETLY True)
endif(VASP_LIBRARY)

if(USE_VASP5)
  set(which_vasp_lib "vasp5")
else(USE_VASP5)
  set(which_vasp_lib "vasp")
endif(USE_VASP5)
FIND_LIBRARY(_VASP_LIBRARY
  ${which_vasp_lib}
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
FIND_PACKAGE_HANDLE_STANDARD_ARGS(escan DEFAULT_MSG _VASP_LIBRARY)
IF(_VASP_LIBRARY)
  set(VASP_LIBRARY ${_VASP_LIBRARY} CACHE PATH "Path to vasp library.")
  SET(VASP_FOUND TRUE)
  unset(_VASP_LIBRARY CACHE)
ENDIF(_VASP_LIBRARY)
