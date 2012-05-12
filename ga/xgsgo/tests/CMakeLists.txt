if(NOT compile_debug_tests)
  add_test(test python ${CMAKE_CURRENT_SOURCE_DIR}/test.py)
endif(NOT compile_debug_tests)

