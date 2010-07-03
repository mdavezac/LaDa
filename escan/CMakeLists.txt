
include_directories(MPI_INCLUDE_PATH)

python_add_module (escan escan.cc python.cc wavefunctions.cc escan_types.f90 fortran.f90 wavefunctions_module.f90)

include_directories(${MPI_INCLUDE_PATH} ${PYTHON_INCLUDE_DIRS})
target_link_libraries(escan ${Boost_PYTHON_LIBRARY} ${MPI_LIBRARIES})
include_directories(${escan_INCLUDE_DIRS})

