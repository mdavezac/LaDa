
include_directories(MPI_INCLUDE_PATH)

python_add_module (_vasp vasp.cc )
include_directories(${MPI_INCLUDE_PATH} ${PYTHON_INCLUDE_DIRS})
target_link_libraries(vasp ${Boost_python_LIBRARY} ${MPI_LIBRARY})

