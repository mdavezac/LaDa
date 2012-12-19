if(NOT compile_debug_tests)

  add_test(gulp_extract_structure python ${CMAKE_CURRENT_SOURCE_DIR}/extract_structure.py)
  add_test(gulp_keywords python ${CMAKE_CURRENT_SOURCE_DIR}/keywords.py)

endif(NOT compile_debug_tests)
