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
              COMPONENTS python mpi mpi_python program_options regex serialization system )
include_directories(${Boost_INCLUDE_DIRS})
#   -- eigen2
find_package( eigen2 REQUIRED )
include_directories(${eigen2_INCLUDE_DIR})
#   -- tinyxml
find_package(tinyxml REQUIRED)
include_directories(${tinyxml_INCLUDE_DIR})


# other packages
find_package(MPI)



# Configure!
if (CMAKE_BUILD_TYPE)
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(LADA_DEBUG TRUE)
  endif (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
endif (CMAKE_BUILD_TYPE)
configure_file (
 "${PROJECT_SOURCE_DIR}/LaDaConfig.h.in"
 "${PROJECT_BINARY_DIR}/LaDaConfig.h"
)


#include root source directory
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(crystal)
add_subdirectory(opt)


option(vasp_wrapper  "Compile vasp wrapper."  off)
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
