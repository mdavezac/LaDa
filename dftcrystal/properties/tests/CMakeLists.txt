if(NOT compile_debug_tests)

  add_test(dftcrystal_properties_keywords python ${CMAKE_CURRENT_SOURCE_DIR}/keywords.py)
  add_test(dftcrystal_properties_bandstructure python ${CMAKE_CURRENT_SOURCE_DIR}/bandstructure.py)

endif(NOT compile_debug_tests)
