set(PYMINIMIZER_PY __init__.py)
set(PYMINIMIZER_HEADERS cgs.hpp  function.hpp  interpolated_gradient.hpp  minimizer.hpp)
set(PYMINIMIZER_SOURCES cgs.cc  interpolated_gradient.cc  minimizer.cc  python.cc)

include_directories(${PYTHON_INCLUDE_DIRS} ${PYUBLAS_INCLUDE_DIRS})
python_add_module (pyminimizer ${PYMINIMIZER_HEADERS} ${PYMINIMIZER_SOURCES})
target_link_libraries(pyminimizer minimizer ${Boost_PYTHON_LIBRARY})
mark_as_advanced(PYTHON_ENABLE_MODULE_pyminimizer)
mark_as_advanced(PYTHON_MODULE_pyenum_BUILD_SHARED)
set_target_properties(pyminimizer PROPERTIES PREFIX "_")
set_target_properties(pyminimizer PROPERTIES OUTPUT_NAME minimizer)

install(TARGETS pyminimizer DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/minimizer)
install(FILES ${PYMINIMIZER_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/minimizer)
