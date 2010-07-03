CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project (LaDa)

# Location where cmake first looks for cmake modules.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

# The version number.
set (LaDa_VERSION_MAJOR 4)
set (LaDa_VERSION_MINOR 5)

enable_language(C)
enable_language(CXX)
enable_language(Fortran)
include(FortranCInterface)
FortranCInterface_HEADER(${PROJECT_BINARY_DIR}/FCMangle.h
                         MACRO_NAMESPACE "FC_"
                         SYMBOL_NAMESPACE "FC_" )


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

# python packages
find_package(PythonInterp)
find_package(PythonLibs)

# other packages
find_package(MPI)



# Configure!
configure_file (
 "${PROJECT_SOURCE_DIR}/LaDaConfig.h.in"
 "${PROJECT_BINARY_DIR}/LaDaConfig.h"
)


#include root source directory
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(crystal)
add_subdirectory(opt)


add_subdirectory(vasp)
