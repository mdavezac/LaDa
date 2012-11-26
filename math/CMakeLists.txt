include_directories(${PYTHON_INCLUDE_DIR} ${NUMPY_INCLUDES})
add_library (pymath SHARED module.cc noopt.cc)
SET_SOURCE_FILES_PROPERTIES(noopt.cc PROPERTIES COMPILE_FLAGS -O0)
set_target_properties(pymath PROPERTIES OUTPUT_NAME math)
set_target_properties(pymath PROPERTIES PREFIX "")
install(TARGETS pymath DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)

