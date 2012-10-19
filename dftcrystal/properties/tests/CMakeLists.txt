if(NOT compile_debug_tests)

  add_test(dftcrystal_properties_keywords python ${CMAKE_CURRENT_SOURCE_DIR}/properties.py)

endif(NOT compile_debug_tests)
