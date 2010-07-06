set(PYENUM_HEADERS bitset.hpp find_all_cells.hpp label_exchange.hpp transform.hpp translation.hpp)
set(PYENUM_SOURCES bitset.cc find_all_cells.cc label_exchange.cc python.cc transform.cc translation.cc)


include_directories(${PYTHON_INCLUDE_DIRS} ${PYUBLAS_INCLUDE_DIRS})
python_add_module (pyenum ${PYENUM_HEADERS} ${PYENUM_SOURCES})
target_link_libraries(pyenum enum ${Boost_PYTHON_LIBRARY})
mark_as_advanced(PYTHON_ENABLE_MODULE_pyenum)
mark_as_advanced(PYTHON_MODULE_pyenum_BUILD_SHARED)
set_target_properties(pyenum PROPERTIES PREFIX "")
set_target_properties(pyenum PROPERTIES OUTPUT_NAME enumeration)

install(TARGETS pyenum DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada)
