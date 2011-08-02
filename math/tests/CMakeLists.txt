if(NOT compile_debug_tests)

  add_executable(symops symops.cc)
  target_link_libraries(symops math)
  ADD_TEST(symops ${CMAKE_CURRENT_BINARY_DIR}/symops)

endif(NOT compile_debug_tests)

