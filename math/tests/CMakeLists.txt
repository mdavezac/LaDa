if(NOT compile_debug_tests)

  set(TYPES 0 1)
  foreach(TYPE ${TYPES})
    add_executable(fuzzy${TYPE} fuzzy.cc)
    set_property( TARGET fuzzy${TYPE} PROPERTY COMPILE_DEFINITIONS
                  LADA_TEST_INCTYPE=${TYPE})
    ADD_TEST(math_fuzzy${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/fuzzy${TYPE})
  endforeach(TYPE)

  set(TYPES 0 1 2 3)
  foreach(TYPE ${TYPES})
    add_executable(fuzzyeigen${TYPE} fuzzyeigen.cc)
    set_property( TARGET fuzzyeigen${TYPE} PROPERTY COMPILE_DEFINITIONS
                  LADA_TEST_INCTYPE=${TYPE})
    ADD_TEST(math_fuzzyeigen${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/fuzzyeigen${TYPE})
  endforeach(TYPE)

  add_executable(symtests symtests.cc)
  ADD_TEST( math_symtests ${CMAKE_CURRENT_BINARY_DIR}/symtests)

  add_library (gruber SHARED "gruber.cc")
  set_target_properties(gruber PROPERTIES PREFIX "_")
  if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    set_property(TARGET gruber PROPERTY COMPILE_DEFINITIONS LADA_LIM=2)
  else("${CMAKE_BUILD_TYPE}" STREQUAL "Debug") 
    set_property(TARGET gruber PROPERTY COMPILE_DEFINITIONS LADA_LIM=5)
  endif("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")

  add_library (smith SHARED "smith_form.cc")
  set_target_properties(smith PROPERTIES PREFIX "_")
  set_target_properties(smith PROPERTIES OUTPUT_NAME "smith")

  add_test(math_basics python ${CMAKE_CURRENT_SOURCE_DIR}/basics.py)
  add_test(math_pygruber python ${CMAKE_CURRENT_SOURCE_DIR}/gruber.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(math_pysmith python ${CMAKE_CURRENT_SOURCE_DIR}/smith.py ${CMAKE_CURRENT_BINARY_DIR})

endif(NOT compile_debug_tests)
