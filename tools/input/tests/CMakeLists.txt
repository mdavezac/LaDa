if(NOT compile_debug_tests)

  add_test(functools_input_tree python ${CMAKE_CURRENT_SOURCE_DIR}/tree.py)
  add_test(functools_input_listblock python ${CMAKE_CURRENT_SOURCE_DIR}/listblock.py)

endif(NOT compile_debug_tests)
