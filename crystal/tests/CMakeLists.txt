if(NOT compile_debug_tests)
endif(NOT compile_debug_tests)

  add_test(atom python ${CMAKE_CURRENT_SOURCE_DIR}/atom.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(atom_memory python ${CMAKE_CURRENT_SOURCE_DIR}/atom_memory.py ${CMAKE_CURRENT_BINARY_DIR})

  add_library (atom_self SHARED "atom_self.cc")
  target_link_libraries(atom_self ${PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARY} crystal mathlib)
  set_target_properties(atom_self PROPERTIES PREFIX "_")
  add_test(atom_self python ${CMAKE_CURRENT_SOURCE_DIR}/atom_self.py ${CMAKE_CURRENT_BINARY_DIR})


  add_test(structure python ${CMAKE_CURRENT_SOURCE_DIR}/structure.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(structure_memory python ${CMAKE_CURRENT_SOURCE_DIR}/structure_memory.py ${CMAKE_CURRENT_BINARY_DIR})

  add_test(utilities python ${CMAKE_CURRENT_SOURCE_DIR}/utilities.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(supercell python ${CMAKE_CURRENT_SOURCE_DIR}/structure_memory.py ${CMAKE_CURRENT_BINARY_DIR})

  #   python_add_module (atom_self "" "atom_self.cc")
  #   set_property(TARGET atom_self PROPERTY COMPILE_DEFINITIONS)
  #   mark_as_advanced(PYTHON_ENABLE_MODULE_atom_self)
  #   mark_as_advanced(PYTHON_MODULE_atom_self_BUILD_SHARED)
  #   target_link_libraries(atom_self pycrystal)
  #   install(TARGETS atom_self DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal)

# MACRO(SET_TYPE increment typename typedef)
#   if(${increment} EQUAL 0)
#     set(${typename} Str)
#     set(${typedef} std::string)
#   elseif(${increment} EQUAL 1)
#     set(${typename} Vec)
#     set(${typedef} std::vector<std::string>)
#   elseif(${increment} EQUAL 2)
#     set(${typename} Set)
#     set(${typedef} std::set<std::string>)
#   endif(${increment} EQUAL 0)
# endmacro(SET_TYPE)

# if(NOT compile_debug_tests)

#   add_executable(utilities utilities.cc)
#   target_link_libraries(utilities crystal)
#   ADD_TEST(utilities ${CMAKE_CURRENT_BINARY_DIR}/utilities)


#   if(loadnsave)
#     set(TYPES  0 1 2)
#     foreach(TYPE ${TYPES})
#       add_executable(lns_atom${TYPE} lns_atom.cc)
#       set_property( TARGET lns_atom${TYPE} PROPERTY COMPILE_DEFINITIONS
#                     LADA_TEST_INCTYPE=${TYPE})
#       target_link_libraries(lns_atom${TYPE} lns)
#       ADD_TEST(lns_atom${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/lns_atom${TYPE})
#     endforeach(TYPE)

#     set(TYPES 0 1 2)
#     set(STRS 1)
#     foreach(STR ${STRS})
#       foreach(TYPE ${TYPES})
#         add_executable(lns_structure${STR}${TYPE} lns_structure.cc)
#         set_property( TARGET lns_structure${STR}${TYPE} PROPERTY COMPILE_DEFINITIONS
#                       LADA_TEST_STRUCTURE=${STR} LADA_TEST_INCTYPE=${TYPE})
#         target_link_libraries(lns_structure${STR}${TYPE} lns)
#         ADD_TEST(lns_structure${STR}${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/lns_structure${STR}${TYPE})
#       endforeach(TYPE)
#     endforeach(STR)
#   endif(loadnsave)

#   add_executable(cellinvariants cellinvariants.cc)
#   target_link_libraries(cellinvariants crystal)
#   ADD_TEST(cellinvariants ${CMAKE_CURRENT_BINARY_DIR}/cellinvariants)

#   add_executable(is_container is_container.cc)
#   target_link_libraries(is_container crystal)
#   ADD_TEST(is_container ${CMAKE_CURRENT_BINARY_DIR}/is_container)

#   set(TYPES  0 1 2 3 4 5 6 7 8)
#   foreach(TYPE ${TYPES})
#     add_executable(init${TYPE} initialize.cc)
#     set_property( TARGET init${TYPE} PROPERTY COMPILE_DEFINITIONS LADA_TEST_INCTYPE=${TYPE})
#     ADD_TEST(init${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/init${TYPE})
#   endforeach(TYPE)

#   set(TYPES 0 1 2 3 4 5 6)
#   foreach(TYPE ${TYPES})
#     add_executable(compare${TYPE} compare.cc)
#     target_link_libraries(compare${TYPE} crystal)
#     set_property( TARGET compare${TYPE} PROPERTY COMPILE_DEFINITIONS
#                   LADA_TEST_INCTYPE=${TYPE})
#     ADD_TEST(compare${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/compare${TYPE})
#   endforeach(TYPE)

#   foreach(INCREMENT 0 1 2)
#     SET_TYPE(INCREMENT TYPENAME TYPE)

#     add_executable(copy_atom${TYPENAME} copy_atom.cc)
#     set_property( TARGET copy_atom${TYPENAME} PROPERTY COMPILE_DEFINITIONS
#       LADA_INCREMENT=${INCREMENT} LADA_TYPE=${TYPE})
#     target_link_libraries(copy_atom${TYPENAME} crystal)
#     ADD_TEST(copy_atom${TYPENAME} ${CMAKE_CURRENT_BINARY_DIR}/copy_atom${TYPENAME})

#     add_executable(copy_structure${TYPENAME} copy_structure.cc)
#     set_property( TARGET copy_structure${TYPENAME} PROPERTY COMPILE_DEFINITIONS
#                   LADA_INCREMENT=${INCREMENT} LADA_TYPE=${TYPE})
#     target_link_libraries(copy_structure${TYPENAME} crystal)
#     ADD_TEST(copy_structure${TYPENAME} ${CMAKE_CURRENT_BINARY_DIR}/copy_structure${TYPENAME})

#     add_executable(supercell_${TYPENAME} supercell.cc)
#     set_property( TARGET supercell_${TYPENAME} PROPERTY COMPILE_DEFINITIONS
#                   LADA_INCREMENT=${INCREMENT} LADA_TYPE=${TYPE})
#     target_link_libraries(supercell_${TYPENAME} crystal)
#     ADD_TEST(supercell_${TYPENAME} ${CMAKE_CURRENT_BINARY_DIR}/supercell_${TYPENAME})

#     add_executable(primitive_${TYPENAME} primitive.cc)
#     target_link_libraries(primitive_${TYPENAME} crystal)
#     if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
#       set_property( TARGET primitive_${TYPENAME} PROPERTY COMPILE_DEFINITIONS
#                     LADA_LIM=5 LADA_INCREMENT=${INCREMENT} LADA_TYPE=${TYPE})
#     else(${CMAKE_BUILD_TYPE} STREQUAL "Debug") 
#       set_property( TARGET primitive_${TYPENAME} PROPERTY COMPILE_DEFINITIONS
#                     LADA_LIM=9 LADA_INCREMENT=${INCREMENT} LADA_TYPE=${TYPE})
#     endif(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
#     ADD_TEST(primitive_${TYPENAME} ${CMAKE_CURRENT_BINARY_DIR}/primitive_${TYPENAME})

#     add_executable(neighbors${TYPENAME} neighbors.cc)
#     set_property( TARGET neighbors${TYPENAME} PROPERTY COMPILE_DEFINITIONS
#                   LADA_INCREMENT=${INCREMENT} LADA_TYPE=${TYPE})
#     target_link_libraries(neighbors${TYPENAME} crystal)
#     ADD_TEST(neighbors${TYPENAME} ${CMAKE_CURRENT_BINARY_DIR}/neighbors${TYPENAME})

#     add_executable(equivalent_structures${TYPENAME} equivalent_structures.cc)
#     set_property( TARGET equivalent_structures${TYPENAME} PROPERTY COMPILE_DEFINITIONS
#                   LADA_INCREMENT=${INCREMENT} LADA_TYPE=${TYPE})
#     target_link_libraries(equivalent_structures${TYPENAME} crystal)
#     ADD_TEST(equivalent_structures${TYPENAME} ${CMAKE_CURRENT_BINARY_DIR}/equivalent_structures${TYPENAME})

#     add_executable(space_group${TYPENAME} space_group.cc)
#     set_property( TARGET space_group${TYPENAME} PROPERTY COMPILE_DEFINITIONS
#                   LADA_INCREMENT=${INCREMENT} LADA_TYPE=${TYPE})
#     target_link_libraries(space_group${TYPENAME} crystal)
#     ADD_TEST(space_group${TYPENAME} ${CMAKE_CURRENT_BINARY_DIR}/space_group${TYPENAME})

#     add_executable(periodic_dnc${TYPENAME} periodic_dnc.cc)
#     set_property( TARGET periodic_dnc${TYPENAME} PROPERTY COMPILE_DEFINITIONS
#                   LADA_INCREMENT=${INCREMENT} LADA_TYPE=${TYPE})
#     target_link_libraries(periodic_dnc${TYPENAME} crystal)
#     ADD_TEST(periodic_dnc${TYPENAME} ${CMAKE_CURRENT_BINARY_DIR}/periodic_dnc${TYPENAME})

#     add_executable(map_sites${TYPENAME} map_sites.cc)
#     set_property( TARGET map_sites${TYPENAME} PROPERTY COMPILE_DEFINITIONS
#                   LADA_INCREMENT=${INCREMENT} LADA_TYPE=${TYPE})
#     target_link_libraries(map_sites${TYPENAME} crystal)
#     ADD_TEST(map_sites${TYPENAME} ${CMAKE_CURRENT_BINARY_DIR}/map_sites${TYPENAME})


#   endforeach(INCREMENT)

# endif(NOT compile_debug_tests)
