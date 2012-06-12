if(NOT compile_debug_tests)

  add_test(functools_makeclass python ${CMAKE_CURRENT_SOURCE_DIR}/makeclass.py ${CMAKE_CURRENT_SOURCE_DIR})

endif(NOT compile_debug_tests)
