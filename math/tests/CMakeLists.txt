if(NOT compile_debug_tests)

  set(TYPES 0 1)
  foreach(TYPE ${TYPES})
    add_executable(fuzzy${TYPE} fuzzy.cc)
    target_link_libraries(fuzzy${TYPE} mathlib)
    set_property( TARGET fuzzy${TYPE} PROPERTY COMPILE_DEFINITIONS
                  LADA_TEST_INCTYPE=${TYPE})
    ADD_TEST(math_fuzzy${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/fuzzy${TYPE})
  endforeach(TYPE)

  set(TYPES 0 1 2 3)
  foreach(TYPE ${TYPES})
    add_executable(fuzzyeigen${TYPE} fuzzyeigen.cc)
    target_link_libraries(fuzzyeigen${TYPE} mathlib)
    set_property( TARGET fuzzyeigen${TYPE} PROPERTY COMPILE_DEFINITIONS
                  LADA_TEST_INCTYPE=${TYPE})
    ADD_TEST(math_fuzzyeigen${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/fuzzyeigen${TYPE})
  endforeach(TYPE)

  add_executable(symtests symtests.cc)
  target_link_libraries(symtests mathlib)
  ADD_TEST( math_symtests ${CMAKE_CURRENT_BINARY_DIR}/symtests)

  add_executable(gruber gruber.cc)
  target_link_libraries(gruber mathlib)
  if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    set_property(TARGET gruber PROPERTY COMPILE_DEFINITIONS LADA_LIM=2)
  else("${CMAKE_BUILD_TYPE}" STREQUAL "Debug") 
    set_property(TARGET gruber PROPERTY COMPILE_DEFINITIONS LADA_LIM=5)
  endif("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  ADD_TEST( math_gruber ${CMAKE_CURRENT_BINARY_DIR}/gruber)

  add_test(math_basics python ${CMAKE_CURRENT_SOURCE_DIR}/basics.py)
  add_test(math_pygruber python ${CMAKE_CURRENT_SOURCE_DIR}/gruber.py)
  add_test(math_pysmith python ${CMAKE_CURRENT_SOURCE_DIR}/smith.py)

  add_executable(smith_form smith_form.cc)
  target_link_libraries(smith_form mathlib)
  ADD_TEST( math_smith_form ${CMAKE_CURRENT_BINARY_DIR}/smith_form)

endif(NOT compile_debug_tests)
