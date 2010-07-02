cmake_minimum_required (2.8)
project (LaDa)

# Location where cmake first looks for cmake modules.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

# The version number.
set (LaDa_VERSION_MAJOR 4)
set (LaDa_VERSION_MINOR 5)

enable_language(CXX)
enable_language(Fortran)

# packages
find_package(stdlibs REQUIRED)
find_package(PythonInterp)
find_package(PythonLibs)


#
configure_file (
 "${PROJECT_SOURCE_DIR}/LaDaConfig.h.in"
 "${PROJECT_BINARY_DIR}/LaDaConfig.h"
)


#include root source directory
include_directories(${LADA_SOURCE_DIR}/w01-intro)

add_subdirectory(crystal)
