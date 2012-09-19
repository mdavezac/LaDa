if(NOT compile_debug_tests)

  add_test(ce_cluster python ${CMAKE_CURRENT_SOURCE_DIR}/cluster.py ${CMAKE_CURRENT_BINARY_DIR})
# add_test(ce_factory python ${CMAKE_CURRENT_SOURCE_DIR}/factory.py ${CMAKE_CURRENT_BINARY_DIR})
  
endif(NOT compile_debug_tests)
