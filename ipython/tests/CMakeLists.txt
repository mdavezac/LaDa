if(NOT compile_debug_tests)

  add_test(explore ipython ${CMAKE_CURRENT_SOURCE_DIR}/explore.py ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(goto ipython ${CMAKE_CURRENT_SOURCE_DIR}/goto.py ${CMAKE_CURRENT_SOURCE_DIR})

endif(NOT compile_debug_tests)
