set(PYVASP_PY methods.py files.py __init__.py kpoints.py launch.py _nosetests.py specie.py)

include_directories(MPI_INCLUDE_PATH)

python_add_module (pyvasp vasp.cc )
include_directories(${MPI_INCLUDE_PATH} ${PYTHON_INCLUDE_DIRS})
target_link_libraries(pyvasp ${Boost_PYTHON_LIBRARY} ${MPI_LIBRARIES})
mark_as_advanced(PYTHON_ENABLE_MODULE_pyvasp)
mark_as_advanced(PYTHON_MODULE_pyvasp_BUILD_SHARED)
set_target_properties(pyvasp PROPERTIES PREFIX "")
set_target_properties(pyvasp PROPERTIES OUTPUT_NAME enumeration)

install(TARGETS pyvasp DESTINATION ${CMAKE_PYINSTALL_PREFIX}/vasp)
install(FILES ${PYVASP_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/vasp)

add_subdirectory(incar)
add_subdirectory(extract)
