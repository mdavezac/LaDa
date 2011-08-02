if(NOT compile_debug_tests)

  set(TYPES 0 1)
  foreach(TYPE ${TYPES})
    add_executable(fuzzy${TYPE} fuzzy.cc)
    target_link_libraries(fuzzy${TYPE} mathlib)
    set_property( TARGET fuzzy${TYPE} PROPERTY COMPILE_DEFINITIONS
                  LADA_TEST_INCTYPE=${TYPE})
    ADD_TEST(fuzzy${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/fuzzy${TYPE})
  endforeach(TYPE)
  set(TYPES 0 1 2 3)
  foreach(TYPE ${TYPES})
    add_executable(fuzzyeigen${TYPE} fuzzyeigen.cc)
    target_link_libraries(fuzzyeigen${TYPE} mathlib)
    set_property( TARGET fuzzyeigen${TYPE} PROPERTY COMPILE_DEFINITIONS
                  LADA_TEST_INCTYPE=${TYPE})
    ADD_TEST(fuzzyeigen${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/fuzzyeigen${TYPE})
  endforeach(TYPE)

endif(NOT compile_debug_tests)
