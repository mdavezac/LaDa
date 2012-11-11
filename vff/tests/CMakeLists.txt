if(NOT compile_debug_tests)

  add_test(vff_node python ${CMAKE_CURRENT_SOURCE_DIR}/node.py)
  add_test(vff_tree python ${CMAKE_CURRENT_SOURCE_DIR}/tree.py)
  add_test(vff_iterator python ${CMAKE_CURRENT_SOURCE_DIR}/iterator.py)
  add_test(vff_represent python ${CMAKE_CURRENT_SOURCE_DIR}/represent.py)
  add_test(vff_ternary python ${CMAKE_CURRENT_SOURCE_DIR}/ternary.py)

endif(NOT compile_debug_tests)
