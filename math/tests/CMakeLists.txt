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

  set(LADA_XS 2 5)
  set(LADA_YS 3 4)
  set(TYPES 0 1)
  foreach(LADA_X ${LADA_XS})
    foreach(LADA_Y ${LADA_YS})
      foreach(TYPE ${TYPES})
        add_executable(serialize${TYPE}-${LADA_X}${LADA_Y} serialize.cc)
        target_link_libraries(serialize${TYPE}-${LADA_X}${LADA_Y} mathlib)
        set_property( TARGET serialize${TYPE}-${LADA_X}${LADA_Y} PROPERTY COMPILE_DEFINITIONS
                      LADA_TEST_INCTYPE=${TYPE} LADA_X=${LADA_X} LADA_Y=${LADA_Y})
        ADD_TEST( math_serialize${TYPE}-${LADA_X}${LADA_Y}
                  ${CMAKE_CURRENT_BINARY_DIR}/serialize${TYPE}-${LADA_X}${LADA_Y})
      endforeach(TYPE)
    endforeach(LADA_Y ${LADA_YS})
  endforeach(LADA_X ${LADA_XS})

  set(TYPES 6 7 8 9)
  foreach(TYPE ${TYPES})
    add_executable(serialize${TYPE} serialize.cc)
    target_link_libraries(serialize${TYPE} mathlib)
    set_property( TARGET serialize${TYPE} PROPERTY COMPILE_DEFINITIONS
                  LADA_TEST_INCTYPE=${TYPE} )
    ADD_TEST( math_serialize${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/serialize${TYPE})
  endforeach(TYPE)

  add_executable(serialize_affine serialize_affine.cc)
  target_link_libraries(serialize_affine mathlib)
  ADD_TEST( math_serialize_affine ${CMAKE_CURRENT_BINARY_DIR}/serialize_affine)

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

  add_executable(smith_form smith_form.cc)
  target_link_libraries(smith_form mathlib)
  ADD_TEST( math_smith_form ${CMAKE_CURRENT_BINARY_DIR}/smith_form)

# add_executable(smith smith.cc)
# target_link_libraries(smith mathlib)
# ADD_TEST( smith ${CMAKE_CURRENT_BINARY_DIR}/smith)

endif(NOT compile_debug_tests)

