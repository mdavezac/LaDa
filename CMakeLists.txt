CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project (LaDa)
set (LaDa_VERSION_MAJOR 4)
set (LaDa_VERSION_MINOR 5)

# Project options
option(disable_executables  "Disable all executables."  off)
mark_as_advanced(disable_executables)

option(pcm_functional      "Include Point Charge Ion functional."       off)
option(ce_functional       "Include Cluster Expansion functional."      off)
option(vasp_functional     "Compile vasp functional wrapper."           off)
option(escan_functional    "Compile escan functional wrapper."          off)
option(vff_functional      "Compile Valence Force Field method."        off)
option(fxseps_functional   "Compile Fixed Lattice Sum of Separables."   off)
option(atompot_functional  "Compile Fixed Lattice Sum of Separables."    on)
option(enum_package        "Compile Lattice Decoration Enumeration."    off)
option(python_bindings     "Include python bindings."                   off)


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
# find_package(PythonInterp)
find_package(PythonLibs)


# required packages
#   -- boost packages
set(Boost_USE_STATIC_LIBS   OFF)
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

if(fxseps_functional)
  if(NOT disable_executables) 
    set(minimizer_libs TRUE)
  endif(NOT disable_executables) 
  set(ce_functional TRUE)
  add_subdirectory(separables)
endif(fxseps_functional)

if(atompot_functional)
  set(python_bindings TRUE)
  add_subdirectory(atomic_potentials)
endif(atompot_functional)

if(minimizer_libs)
  find_package(gsl)
  set(LADA_WITH_GSL GSL_FOUND)
  find_package(minuit2)
  set(LADA_WITH_MINUIT2 MINUIT2_FOUND)
  add_subdirectory(minimizer)
endif(minimizer_libs)

if(python_bindings)
  if(PYTHONLIBS_FOUND)
    add_subdirectory(python)
  else(PYTHONLIBS_FOUND)
    message(SEND_ERROR "Python libraries not found."
                       "Please set PYTHON_INCLUDE_DIR and PYTHON_LIBRARY" )
    set(DOPYTHON False)
  endif(PYTHONLIBS_FOUND)
endif(python_bindings)

# Configure file.
if (CMAKE_BUILD_TYPE)
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(LADA_DEBUG TRUE)
  endif (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
endif (CMAKE_BUILD_TYPE)
configure_file (
 "${PROJECT_SOURCE_DIR}/LaDaConfig.h.in"
 "${PROJECT_BINARY_DIR}/LaDaConfig.h"
)

