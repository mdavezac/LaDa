if(NOT compile_debug_tests)

  add_executable(utilities utilities.cc)
  target_link_libraries(utilities crystal)
  ADD_TEST(utilities ${CMAKE_CURRENT_BINARY_DIR}/utilities)


  if(loadnsave)
    set(TYPES  0 1 2)
    foreach(TYPE ${TYPES})
      add_executable(lns_atom${TYPE} lns_atom.cc)
      set_property( TARGET lns_atom${TYPE} PROPERTY COMPILE_DEFINITIONS
                    LADA_TEST_INCTYPE=${TYPE})
      target_link_libraries(lns_atom${TYPE} lns)
      ADD_TEST(lns_atom${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/lns_atom${TYPE})
    endforeach(TYPE)

    set(TYPES 0 1 2)
    set(STRS 1)
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

  add_executable(cellinvariants cellinvariants.cc)
  target_link_libraries(cellinvariants crystal)
  ADD_TEST(cellinvariants ${CMAKE_CURRENT_BINARY_DIR}/cellinvariants)

  add_executable(is_container is_container.cc)
  target_link_libraries(is_container crystal)
  ADD_TEST(is_container ${CMAKE_CURRENT_BINARY_DIR}/is_container)

  set(TYPES  0 1 2 3 4 5 6 7 8)
  foreach(TYPE ${TYPES})
    add_executable(init${TYPE} initialize.cc)
    set_property( TARGET init${TYPE} PROPERTY COMPILE_DEFINITIONS LADA_TEST_INCTYPE=${TYPE})
    ADD_TEST(init${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/init${TYPE})
  endforeach(TYPE)

  set(TYPES 0 1 2 3 4 5 6)
  foreach(TYPE ${TYPES})
    add_executable(compare${TYPE} compare.cc)
    target_link_libraries(compare${TYPE} crystal)
    set_property( TARGET compare${TYPE} PROPERTY COMPILE_DEFINITIONS
                  LADA_TEST_INCTYPE=${TYPE})
    ADD_TEST(compare${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/compare${TYPE})
  endforeach(TYPE)

  set(TYPES "set" "vector")
  foreach(TYPE ${TYPES})

    add_executable(supercell_${TYPE} supercell.cc)
    set_property(TARGET supercell_${TYPE} PROPERTY COMPILE_DEFINITIONS LADA_TYPE=std::${TYPE}<std::string>)
    target_link_libraries(supercell_${TYPE} crystal)
    ADD_TEST(supercell_${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/supercell_${TYPE})

    add_executable(primitive_${TYPE} primitive.cc)
    target_link_libraries(primitive_${TYPE} crystal)
    if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
      set_property( TARGET primitive_${TYPE} PROPERTY COMPILE_DEFINITIONS
                    LADA_LIM=5 LADA_TYPE=std::${TYPE}<std::string> )
    else(${CMAKE_BUILD_TYPE} STREQUAL "Debug") 
      set_property( TARGET primitive_${TYPE} PROPERTY COMPILE_DEFINITIONS
                    LADA_LIM=9 LADA_TYPE=std::${TYPE}<std::string> )
    endif(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    ADD_TEST(primitive_${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/primitive_${TYPE})

    add_executable(neighbors_${TYPE} neighbors.cc)
    set_property(TARGET neighbors_${TYPE} PROPERTY COMPILE_DEFINITIONS LADA_TYPE=std::${TYPE}<std::string>)
    target_link_libraries(neighbors_${TYPE} crystal)
    ADD_TEST(supercell ${CMAKE_CURRENT_BINARY_DIR}/neighbors_${TYPE})

    add_executable(equivalent_structures_${TYPE} equivalent_structures.cc)
    set_property( TARGET equivalent_structures_${TYPE} PROPERTY COMPILE_DEFINITIONS
                  LADA_TYPE=std::${TYPE}<std::string> )
    target_link_libraries(equivalent_structures_${TYPE} crystal)
    ADD_TEST(equivalent_structures_${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/equivalent_structures_${TYPE})

    add_executable(space_group_${TYPE} space_group.cc)
    set_property(TARGET space_group_${TYPE} PROPERTY COMPILE_DEFINITIONS LADA_TYPE=std::${TYPE}<std::string>)
    target_link_libraries(space_group_${TYPE} crystal)
    ADD_TEST(space_group_${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/space_group_${TYPE})

  endforeach(TYPE ${TYPES})

endif(NOT compile_debug_tests)

    set(TYPE "set")
    add_executable(map_sites_${TYPE} map_sites.cc)
    set_property(TARGET map_sites_${TYPE} PROPERTY COMPILE_DEFINITIONS LADA_TYPE=std::${TYPE}<std::string>)
    target_link_libraries(map_sites_${TYPE} crystal)
    ADD_TEST(map_sites_${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/map_sites_${TYPE})

