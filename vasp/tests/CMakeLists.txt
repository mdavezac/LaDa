if(NOT compile_debug_tests)
  

endif(NOT compile_debug_tests)
add_test(specie python ${CMAKE_CURRENT_SOURCE_DIR}/species.py ${CMAKE_CURRENT_SOURCE_DIR})

