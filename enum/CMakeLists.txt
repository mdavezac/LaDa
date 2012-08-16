set (PYCE_PY __init__.py transforms.py defects.py)

install(FILES ${PYCE_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/enum)

include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
add_library (pyenum SHARED module.cc ndimiterator.cc ndimiterator.h fciterator.cc fciterator.h)
set_target_properties(pyenum PROPERTIES OUTPUT_NAME cppwrappers)
set_target_properties(pyenum PROPERTIES PREFIX "")
target_link_libraries(pyenum ${PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARY} mathlib)
install(TARGETS pyenum DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/enum)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
