if(NOT compile_debug_tests)
  add_test(test python ${CMAKE_CURRENT_SOURCE_DIR}/test.py
                       ${CMAKE_CURRENT_BINARY_DIR}/../../process/test/pifunc)
endif(NOT compile_debug_tests)

