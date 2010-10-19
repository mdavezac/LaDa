set(PYVASP_PY methods.py files.py __init__.py kpoints.py launch.py _nosetests.py specie.py)

include_directories(MPI_INCLUDE_PATH)

if(VASP_FOUND)
  python_add_module (pyvasp vasp.cc )
  include_directories(${PYTHON_INCLUDE_DIRS})
  target_link_libraries(pyvasp ${VASP_LIBRARY} ${Boost_PYTHON_LIBRARY})
  if(MPI_FOUND)
     include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(pyvasp ${MPI_LIBRARIES} ${Boost_MPI_LIBRARY} )
  endif(MPI_FOUND)
  mark_as_advanced(PYTHON_ENABLE_MODULE_pyvasp)
  mark_as_advanced(PYTHON_MODULE_pyvasp_BUILD_SHARED)
  set_target_properties(pyvasp PROPERTIES PREFIX "_")
  set_target_properties(pyvasp PROPERTIES OUTPUT_NAME vasp)

  install(TARGETS pyvasp DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/vasp)
endif(VASP_FOUND)
install(FILES ${PYVASP_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/vasp)

add_subdirectory(incar)
add_subdirectory(extract)
