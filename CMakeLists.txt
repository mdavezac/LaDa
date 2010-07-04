CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project (LaDa)
set (LaDa_VERSION_MAJOR 4)
set (LaDa_VERSION_MINOR 5)

# Project advanced options 
option(disable_executables  "Disable all executables."  off)
mark_as_advanced(disable_executables)
option(noinstall_headers    "Will not install header files." on)
mark_as_advanced(noinstall_headers)
option(load_python_stuff_globally "Whether python should load escan/vasp globally." off)
mark_as_advanced(load_python_stuff_globally)
option(fucking_cray "Cray compute nodes don't have access to /dev/std*." off)
mark_as_advanced(fucking_cray)
option(use_gsl "Whether to use GSL as minimizers." on)
mark_as_advanced(use_gsl)
option(use_minuit2 "Whether to use Minuit2 as minimizers." on)
mark_as_advanced(use_minuit2)

# Project package options 
option(pcm_functional         "Include Point Charge Ion functional."       off)
option(ce_functional          "Include Cluster Expansion functional."      off)
option(vasp_functional        "Compile vasp functional wrapper."           off)
option(escan_functional       "Compile escan functional wrapper."          off)
option(vff_functional         "Compile Valence Force Field method."        off)
option(separables_functional  "Compile Fixed Lattice Sum of Separables."   off)
option(atompot_functional     "Compile Fixed Lattice Sum of Separables."   off)
option(enum_package           "Compile Lattice Decoration Enumeration."    off)
option(python_bindings        "Include python bindings."                   off)


# Location where cmake first looks for cmake modules.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)


# Languages
enable_language(C)
enable_language(CXX)
enable_language(Fortran)
include(FortranCInterface)
FortranCInterface_HEADER(${PROJECT_BINARY_DIR}/FCMangle.h
                         MACRO_NAMESPACE "FC_"
                         SYMBOL_NAMESPACE "FC_" )
find_package(PythonInterp)
find_package(PythonLibs)
# look for python version.
if(NOT PYTHON_VERSION AND PYTHONINTERP_FOUND)
  execute_process( 
    COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print \"%i.%i\" % sys.version_info[:2]"
    OUTPUT_VARIABLE PYTHON_VERSION
  )
  if( PYTHON_VERSION )
    string (STRIP ${PYTHON_VERSION} PYTHON_VERSION)
    set(PYTHON_VERSION ${PYTHON_VERSION} CACHE STRING "Version of the Python interpreter.")
    mark_as_advanced(PYTHON_VERSION)
    MESSAGE(STATUS "[Python] Version: ${PYTHON_VERSION}")
  else( PYTHON_VERSION )
    MESSAGE(STATUS "Could not determine python version.")
  endif( PYTHON_VERSION )
endif(NOT PYTHON_VERSION AND PYTHONINTERP_FOUND)
# sets python site-packages
if(NOT CMAKE_PYINSTALL_PREFIX)
  execute_process( 
    COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print get_python_lib()"
    OUTPUT_VARIABLE CMAKE_PYINSTALL_PREFIX
  )
  if( CMAKE_PYINSTALL_PREFIX )
    string (STRIP ${CMAKE_PYINSTALL_PREFIX} CMAKE_PYINSTALL_PREFIX)
    set(CMAKE_PYINSTALL_PREFIX ${CMAKE_PYINSTALL_PREFIX} CACHE STRING "Version of the Python interpreter.")
    mark_as_advanced(CMAKE_PYINSTALL_PREFIX)
    MESSAGE(STATUS "[Python] Version: ${CMAKE_PYINSTALL_PREFIX}")
  else( CMAKE_PYINSTALL_PREFIX )
    set( CMAKE_PYINSTALL_PREFIX lib64/python${PYTHON_VERSION}/site-packages/ 
         CACHE PATH "Python modules will be installed here." )
  endif( CMAKE_PYINSTALL_PREFIX )
endif(NOT CMAKE_PYINSTALL_PREFIX)


# required packages
#   -- boost packages
set(BOOST_ROOT   $ENV{HOME}/usr)
find_package( Boost 1.36.0 REQUIRED 
              COMPONENTS python mpi mpi_python program_options regex serialization system filesystem )
include_directories(${Boost_INCLUDE_DIRS})
#   -- eigen2
find_package( eigen2 REQUIRED )
include_directories(${eigen2_INCLUDE_DIR})
#   -- tinyxml
find_package(tinyxml REQUIRED)


# other packages
find_package(MPI)
if( MPI_FOUND )
  set(LADA_MPI TRUE)
  set(_MPI TRUE)
endif( MPI_FOUND )
  

# Look for numpy. Not all packages require NUMPY, however, it is usefull for
# most python bindings.
find_package(numpy)



#include root source directory
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_BINARY_DIR}")

# the most basic stuff to compile/include
add_subdirectory(mpi)
add_subdirectory(opt)
add_subdirectory(math)
add_subdirectory(crystal)


if( vasp_functional )
  set(DOVASP True)
  if( NOT MPI_FOUND )
    message(SEND_ERROR "Vasp interface requires MPI.")
    set(DOVASP False)
  endif( NOT MPI_FOUND )
  if( DOVASP )
    add_subdirectory(vasp)
    set(python_bindings TRUE)
  endif( DOVASP )
endif( vasp_functional )


if( escan_functional )
  set(DOESCAN True)
  if( NOT MPI_FOUND )
    message(SEND_ERROR "Escan interface requires MPI.")
    set(DOESCAN False)
  endif( NOT MPI_FOUND )
  if( NOT NUMPY_FOUND )
    message(SEND_ERROR "Escan interface requires numpy and friends.")
    set(DOESCAN False)
  endif( NOT NUMPY_FOUND )
  find_package(escan REQUIRED)
  if( DOESCAN )
    add_subdirectory(escan)
    set(vff_functional TRUE)
    set(python_bindings TRUE)
  endif( DOESCAN )
endif( escan_functional )

if( vff_functional )
  set(minimizer_libs TRUE)
  add_subdirectory(vff)
endif( vff_functional )

if( ce_functional )
  set(minimizer_libs TRUE)
  add_subdirectory(ce)
endif( ce_functional )

if(pcm_functional)
  set(minimizer_libs TRUE)
  add_subdirectory(pcm)
endif(pcm_functional)

if(enum_package)
  set(python_bindings TRUE)
  add_subdirectory(enum)
endif(enum_package)

if(separables_functional)
  if(NOT disable_executables) 
    set(minimizer_libs TRUE)
  endif(NOT disable_executables) 
  set(ce_functional TRUE)
  add_subdirectory(separables)
endif(separables_functional)

if(atompot_functional)
  set(python_bindings TRUE)
  add_subdirectory(atomic_potentials)
endif(atompot_functional)

if(minimizer_libs)
  if(use_gsl) 
    find_package(gsl)
  endif(use_gsl) 
  set(LADA_WITH_GSL GSL_FOUND)
  if(use_minuit2) 
    find_package(minuit2)
  endif(use_minuit2) 
  set(LADA_WITH_MINUIT2 MINUIT2_FOUND)
  add_subdirectory(minimizer)
endif(minimizer_libs)

if(python_bindings)
  if(PYTHONLIBS_FOUND)
    add_subdirectory(python)
    set(LADA_DO_PYTHON TRUE)
  else(PYTHONLIBS_FOUND)
    message(SEND_ERROR "Python libraries not found."
                       "Please set PYTHON_INCLUDE_DIR and PYTHON_LIBRARY" )
    set(DOPYTHON False)
  endif(PYTHONLIBS_FOUND)
endif(python_bindings)

# Configure file.
if (load_python_stuff_globally)
  set(load_python_stuff_globally true)
else (load_python_stuff_globally)
  set(load_python_stuff_globally false)
endif (load_python_stuff_globally)
if (CMAKE_BUILD_TYPE)
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(LADA_DEBUG TRUE)
  endif (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
endif (CMAKE_BUILD_TYPE)
configure_file (
 "${PROJECT_SOURCE_DIR}/LaDaConfig.h.in"
 "${PROJECT_BINARY_DIR}/LaDaConfig.h"
)


# #  Adds a global library. Also used by python.
# if(python_bindings OR NOT noinstall_headers) 
#   # We want to compile with fPIC in this case.
#   include(CheckCXXCompilerFlag)
#   CHECK_CXX_COMPILER_FLAG(-fPIC CCOMPILER_HAS_FPIC)
#   if(COMPILER_HAS_FPIC)
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
#     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
#   endif(COMPILER_HAS_FPIC)
#   
#   # now creates library.
#   add_library(lada "")
#   target_link_libraries(lada ${Boost_PYTHON_LIBRARY} crystal opt math)
#   if(MPI_FOUND)
#     include_directories(${MPI_INCLUDE_PATH})
#     target_link_libraries(lada mpi)
#   endif(MPI_FOUND)
#   if(vff_functional)
#     target_link_libraries(lada vff)
#   endif(vff_functional)
#   if(ce_functional)
#     target_link_libraries(lada ce)
#   endif(ce_functional)
#   if(atompot_functional)
#     target_link_libraries(lada atompot)
#   endif(atompot_functional)
#   if(pcm_functional)
#     target_link_libraries(lada pcm)
#   endif(pcm_functional)
#   if(enum_package)
#     target_link_libraries(lada enum)
#   endif(enum_package)
#   if(minimizer_libs)
#     target_link_libraries(lada minimizer)
#   endif(minimizer_libs)

#   if(NOT noinstall_headers)
#     install(FILES ${PYLADA_HEADERS} DESTINATION include/LaDa/python)
#   endif(NOT noinstall_headers)
# endif(python_bindings OR NOT noinstall_headers) 
