CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project (LaDa)
set (LaDa_VERSION_MAJOR 5)
set (LaDa_VERSION_MINOR 0)

# Project advanced options 
option(disable_executables  "Disable all executables."  off)
mark_as_advanced(disable_executables)
option(noinstall_headers    "Will not install header files." on)
mark_as_advanced(noinstall_headers)
option(load_python_stuff_globally "Whether python should load escan/vasp globally." off)
mark_as_advanced(load_python_stuff_globally)
option(FUCKING_CRAY "Cray compute nodes don't have access to /dev/std*." off)
mark_as_advanced(FUCKING_CRAY)
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
option(jobs_package           "Installs jobs python module."               off)
option(ga_package             "Installs GA python module."                 off)
option(python_bindings        "Include python bindings."                   off)
option(ipython_interface      "Include IPython interface."                 off)
option(ladabase_package       "Install ladabase package for pymongo."      off)


# Location where cmake first looks for cmake modules.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

# Location where cmake looks for libraries
# adds environment SOFTWARE_ROOT if it exists.
# adds SOFTWARE_ROOT if it exists.
# adds environment LD_LIBRARY_PATH if it exists.
if (NOT "$ENV{SOFTWARE_ROOT}" STREQUAL "")
  set(CMAKE_LIBRARY_PATH $ENV{SOFTWARE_ROOT}/lib $ENV{SOFTWARE_ROOT}/lib64  ${CMAKE_LIBRARY_PATH})
  set(CMAKE_INCLUDE_PATH $ENV{SOFTWARE_ROOT}/include)
else (NOT "$ENV{SOFTWARE_ROOT}" STREQUAL "")
  set(CMAKE_LIBRARY_PATH $ENV{HOME}/usr/lib $ENV{HOME}/usr/lib64  ${CMAKE_LIBRARY_PATH})
  set(CMAKE_INCLUDE_PATH $ENV{HOME}/usr/include)
endif (NOT "$ENV{SOFTWARE_ROOT}" STREQUAL "")
if(SOFTWARE_ROOT)
  set(CMAKE_LIBRARY_PATH ${SOFTWARE_ROOT}/lib ${SOFTWARE_ROOT}/lib64  ${CMAKE_LIBRARY_PATH})
  set(CMAKE_INCLUDE_PATH ${SOFTWARE_ROOT}/include)
endif (SOFTWARE_ROOT)
if (NOT "$ENV{LD_LIBRARY_PATH}" STREQUAL "")
  string(REPLACE ":" ";" _LDLIBPATH $ENV{LD_LIBRARY_PATH})
  set(CMAKE_LIBRARY_PATH ${_LDLIBPATH} ${CMAKE_LIBRARY_PATH})
endif (NOT "$ENV{LD_LIBRARY_PATH}" STREQUAL "")
  


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
  if(PYTHON_VERSION)
    set( CMAKE_PYINSTALL_PREFIX lib64/python${PYTHON_VERSION}/site-packages/ 
         CACHE PATH "Python modules will be installed here." )
  else(PYTHON_VERSION)
    execute_process( 
      COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print get_python_lib()"
      OUTPUT_VARIABLE CMAKE_PYINSTALL_PREFIX
    )
    if( CMAKE_PYINSTALL_PREFIX )
      string (STRIP ${CMAKE_PYINSTALL_PREFIX} CMAKE_PYINSTALL_PREFIX)
      set(CMAKE_PYINSTALL_PREFIX ${CMAKE_PYINSTALL_PREFIX} CACHE STRING "Version of the Python interpreter.")
      mark_as_advanced(CMAKE_PYINSTALL_PREFIX)
      MESSAGE(STATUS "[Python] Version: ${CMAKE_PYINSTALL_PREFIX}")
    endif( CMAKE_PYINSTALL_PREFIX )
  endif(PYTHON_VERSION)
endif(NOT CMAKE_PYINSTALL_PREFIX)


# required packages
#   -- boost packages
if(NOT BOOST_ROOT)
  if($ENV{BOOST_ROOT})
    set(BOOST_ROOT $ENV{BOOST_ROOT} CACHE PATH "Path to boost libraries and headers.")
  else($ENV{BOOST_ROOT})
    set(BOOST_ROOT $ENV{HOME}/usr CACHE PATH "Path to boost libraries and headers.")
  endif($ENV{BOOST_ROOT})
endif(NOT BOOST_ROOT)
find_package( Boost 1.36.0 REQUIRED 
              COMPONENTS python mpi mpi_python program_options regex serialization system filesystem )
include_directories(${Boost_INCLUDE_DIRS})
#   -- eigen2
find_package( eigen2 REQUIRED )
include_directories(${eigen2_INCLUDE_DIR})
#   -- tinyxml
find_package(tinyxml REQUIRED)


# MPI packages. Cray's fucking things up as usual
if("$ENV{NERSC_HOST}" STREQUAL "hopper")
  MESSAGE(STATUS "Working on hopper.")
  MESSAGE(STATUS "Crays being difficult beasts, options are set directly in CMakeLists.txt")
  set(FORTRAN_RTLIBRARIES /opt/pgi/10.3.0/linux86-64/10.3/libso/libpgf90.so      
                          /opt/pgi/10.3.0/linux86-64/10.3/libso/libpgf90_rpm1.so 
                          /opt/pgi/10.3.0/linux86-64/10.3/libso/libpgf902.so     
                          /opt/pgi/10.3.0/linux86-64/10.3/libso/libpgf90rtl.so   
                          /opt/pgi/10.3.0/linux86-64/10.3/libso/libpgftnrtl.so  
                          /opt/pgi/10.3.0/linux86-64/10.3/libso/libpgc.so  )
  set(MPI_INCLUDE_PATH "/opt/cray/mpt/4.1.0.1/xt/seastar/mpich2-pgi64/include/")
  set(MPI_COMPILER "CC")
  set(MPIEXEC "aprun")
  set(MPI_FOUND TRUE)
  set(FUCKING_CRAY TRUE)
  set(load_python_stuff_globally TRUE)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DFUCKING_CRAY")
elseif("$ENV{NERSC_HOST}" STREQUAL "hopper2")
  MESSAGE(STATUS "Working on hopper2.")
  MESSAGE(STATUS "Crays being difficult beasts, options are set directly in CMakeLists.txt")
  set(MPI_INCLUDE_PATH "/tmp")
  set(MPI_COMPILER "CC")
  set(MPIEXEC "aprun")
  set(MPI_FOUND TRUE)
  set(FUCKING_CRAY TRUE)
  set(load_python_stuff_globally FALSE)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DFUCKING_CRAY")
else("$ENV{NERSC_HOST}" STREQUAL "hopper")
  find_package(MPI)
  FortranCInterface_VERIFY()
endif("$ENV{NERSC_HOST}" STREQUAL "hopper")

if( MPI_FOUND )
  set(LADA_MPI TRUE)
endif( MPI_FOUND )
  




#include root source directory
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_BINARY_DIR}")

# the most basic stuff to compile/include
add_subdirectory(opt)
add_subdirectory(math)
add_subdirectory(crystal)
add_subdirectory(physics)

# checks for package dependencies.
if(ladabase_package)
  set(ipython_interface TRUE)
  set(python_bindings TRUE)
  set(vasp_functional TRUE)
endif(ladabase_package)
if(escan_functional)
  set(vff_functional TRUE)
endif(escan_functional)
if(escan_functional OR vff_functional OR pcm_functional)
  set(minimizer_libs TRUE)
endif(escan_functional OR vff_functional OR pcm_functional)
if(escan_functional OR vasp_functional OR enum_package OR atompot_functional OR ipython_interface)
  set(python_bindings TRUE)
endif(escan_functional OR vasp_functional OR enum_package OR atompot_functional OR ipython_interface)
if(ipython_interface)
  set(jobs_package on)
endif(ipython_interface)
if(enum_package)
  set(ce_functional on)
endif(enum_package)
if(NOT disable_executables) 
  if(ce_functional OR separables_functional)
    set(minimizer_libs on)
  endif(ce_functional OR separables_functional)
endif(NOT disable_executables) 

# now includes packages.
if(minimizer_libs)
  if(use_gsl) 
    find_package(gsl)
    set(LADA_WITH_GSL GSL_FOUND)
  else(use_gsl) 
    set(LADA_WITH_GSL FALSE)
  endif(use_gsl) 
  if(use_minuit2) 
    find_package(minuit2)
    set(LADA_WITH_MINUIT2 MINUIT2_FOUND)
  else(use_minuit2) 
    set(LADA_WITH_MINUIT2 FALSE)
  endif(use_minuit2) 
  add_subdirectory(minimizer)
endif(minimizer_libs)

if(python_bindings)
  if(PYTHONLIBS_FOUND)
    add_subdirectory(python)
    add_subdirectory(periodic_table)
    set(LADA_DO_PYTHON TRUE)
  else(PYTHONLIBS_FOUND)
    message(SEND_ERROR "Python libraries not found."
                       "Please set PYTHON_INCLUDE_DIR and PYTHON_LIBRARY" )
    set(DOPYTHON False)
  endif(PYTHONLIBS_FOUND)
  find_package(numpy REQUIRED)
endif(python_bindings)

if( vasp_functional )
  set(DOVASP True)
  add_subdirectory(vasp)
endif( vasp_functional )

if(ladabase_package)
  find_package(pymongo REQUIRED)
  add_subdirectory(ladabase)
endif(ladabase_package)


set(LADA_DO_ESCAN FALSE)
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
  find_package(escan)
  if( DOESCAN )
    add_subdirectory(escan)
    set(LADA_DO_ESCAN TRUE)
  endif( DOESCAN )
endif( escan_functional )

if( vff_functional )
  add_subdirectory(vff)
endif( vff_functional )

if( ce_functional )
  add_subdirectory(ce)
endif( ce_functional )

if(pcm_functional)
  add_subdirectory(pcm)
endif(pcm_functional)

if(enum_package)
  add_subdirectory(enum)
endif(enum_package)

if(separables_functional)
  set(ce_functional on)
  add_subdirectory(separables)
endif(separables_functional)

if(atompot_functional)
  add_subdirectory(atomic_potentials)
endif(atompot_functional)

if(ipython_interface)
  add_subdirectory(ipython)
endif(ipython_interface)

if(jobs_package)
  add_subdirectory(jobs)
endif(jobs_package)
if(ga_package)
  add_subdirectory(ga)
endif(ga_package)

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

