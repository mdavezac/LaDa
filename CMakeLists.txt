CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project (LaDa)

# Location where cmake first looks for cmake modules.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

# The version number.
set (LaDa_VERSION_MAJOR 4)
set (LaDa_VERSION_MINOR 5)

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
  




#include root source directory
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(mpi)
add_subdirectory(opt)
add_subdirectory(math)
add_subdirectory(crystal)

option(disable_executables  "Disable all executables."  off)
mark_as_advanced(disable_executables)
option(vasp_wrapper  "Compile vasp wrapper."  off)
option(escan_functional  "Compile escan python wrapper."  on)
option(vff_functional  "Compile Valence Force Field method."  off)

if( vasp_wrapper )
  set(DOVASP True)
  if( NOT MPI_FOUND )
    message(SEND_ERROR "Vasp interface requires MPI.")
    set(DOVASP False)
  endif( NOT MPI_FOUND )
  if( NOT PYTHONLIBS_FOUND )
    message(SEND_ERROR "Vasp interface requires Python libraries."
                       "Please set PYTHON_INCLUDE_DIR and PYTHON_LIBRARY" )
    set(DOVASP False)
  endif( NOT PYTHONLIBS_FOUND )
  if( DOVASP )
    add_subdirectory(vasp)
  endif( DOVASP )
endif( vasp_wrapper )


if( escan_functional )
  set(DOESCAN True)
  if( NOT MPI_FOUND )
    message(SEND_ERROR "Vasp interface requires MPI.")
    set(DOESCAN False)
  endif( NOT MPI_FOUND )
  if( NOT PYTHONLIBS_FOUND )
    message(SEND_ERROR "Vasp interface requires Python libraries."
                       "Please set PYTHON_INCLUDE_DIR and PYTHON_LIBRARY" )
    set(DOESCAN False)
  endif( NOT PYTHONLIBS_FOUND )
  find_package(escan REQUIRED)
  find_package(numpy REQUIRED)
  if( DOESCAN )
    add_subdirectory(escan)
    set(vff_functional TRUE)
  endif( DOESCAN )
endif( escan_functional )

if( vff_functional )
  set(minimizer_libs TRUE)
  add_subdirectory(vff)
endif( vff_functional )

if(minimizer_libs)
  find_package(gsl)
  find_package(minuit2)
  add_subdirectory(minimizer)
endif(minimizer_libs)

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

