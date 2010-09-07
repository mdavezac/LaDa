set(PYPHYSICS_SOURCES python.cc)


include_directories(${PYTHON_INCLUDE_DIRS})
python_add_module (pyphysics ${PYPHYSICS_HEADERS} ${PYPHYSICS_SOURCES})
target_link_libraries(pyphysics ${Boost_PYTHON_LIBRARY})
mark_as_advanced(PYTHON_ENABLE_MODULE_pyphysics)
mark_as_advanced(PYTHON_MODULE_pyphysics_BUILD_SHARED)
set_target_properties(pyphysics PROPERTIES PREFIX "_")
set_target_properties(pyphysics PROPERTIES OUTPUT_NAME physics)
install(TARGETS pyphysics DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/physics)
install(FILES __init__.py DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/physics)

if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
endif(MPI_FOUND)
