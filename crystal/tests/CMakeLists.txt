if(NOT compile_debug_tests)

  set(TYPES  0 1 2 3)
  set(STRS StructureData TemplateStructure)
  foreach(STR ${STRS})
    foreach(TYPE ${TYPES})
      add_executable(init${STR}${TYPE} initialize.cc)
      set_property( TARGET init${STR}${TYPE} PROPERTY COMPILE_DEFINITIONS
                    LADA_TEST_STRUCTURE=${STR} LADA_TEST_INCTYPE=${TYPE})
      ADD_TEST(init${STR}${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/init${STR}${TYPE})
    endforeach(TYPE)
  endforeach(STR)


  add_executable(utilities utilities.cc)
  target_link_libraries(utilities crystal)
  ADD_TEST(utilities ${CMAKE_CURRENT_BINARY_DIR}/utilities)

  set(TYPES 0 1 2 3 4 5)
  foreach(TYPE ${TYPES})
    add_executable(compare${TYPE} compare.cc)
    target_link_libraries(compare${TYPE} crystal)
    set_property( TARGET compare${TYPE} PROPERTY COMPILE_DEFINITIONS
                  LADA_TEST_INCTYPE=${TYPE})
    ADD_TEST(compare${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/compare${TYPE})
  endforeach(TYPE)

  set(TYPES 0 1 2 3 4 5)
  foreach(TYPE ${TYPES})
    add_executable(is_container${TYPE} is_container.cc)
    target_link_libraries(is_container${TYPE} crystal)
    set_property( TARGET is_container${TYPE} PROPERTY COMPILE_DEFINITIONS
                  LADA_TEST_INCTYPE=${TYPE})
    ADD_TEST(is_container${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/is_container${TYPE})
  endforeach(TYPE)

  if(loadnsave)
    set(TYPES  0 1)
    foreach(TYPE ${TYPES})
      add_executable(lns_atom${TYPE} lns_atom.cc)
      set_property( TARGET lns_atom${TYPE} PROPERTY COMPILE_DEFINITIONS
                    LADA_TEST_INCTYPE=${TYPE})
      target_link_libraries(lns_atom${TYPE} lns)
      ADD_TEST(lns_atom${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/lns_atom${TYPE})
    endforeach(TYPE)

    set(TYPES  0 1)
    set(STRS 0 1)
    foreach(STR ${STRS})
      foreach(TYPE ${TYPES})
        add_executable(lns_structure${STR}${TYPE} lns_structure.cc)
        set_property( TARGET lns_structure${STR}${TYPE} PROPERTY COMPILE_DEFINITIONS
                      LADA_TEST_STRUCTURE=${STR} LADA_TEST_INCTYPE=${TYPE})
        target_link_libraries(lns_structure${STR}${TYPE} lns)
        ADD_TEST(lns_structure${STR}${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/lns_structure${STR}${TYPE})
      endforeach(TYPE)
    endforeach(STR)
  endif(loadnsave)

  add_executable(supercell supercell.cc)
  target_link_libraries(supercell crystal)
  ADD_TEST(supercell ${CMAKE_CURRENT_BINARY_DIR}/supercell)

endif(NOT compile_debug_tests)

  add_executable(equivalent_structures equivalent_structures.cc)
  target_link_libraries(equivalent_structures crystal)
  ADD_TEST(equivalent_structures ${CMAKE_CURRENT_BINARY_DIR}/equivalent_structures)

# add_executable(primitive primitive.cc)
# target_link_libraries(primitive crystal)
# ADD_TEST(primitive ${CMAKE_CURRENT_BINARY_DIR}/primitive)

