  add_executable (pifunc functional.cc)
if(NOT compile_debug_tests)
  # add_test(test python ${CMAKE_CURRENT_SOURCE_DIR}/test.py ${CMAKE_CURRENT_BINARY_DIR}/pifunc)
endif(NOT compile_debug_tests)

