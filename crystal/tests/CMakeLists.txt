if(NOT compile_debug_tests)

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
  add_test(primitive python ${CMAKE_CURRENT_SOURCE_DIR}/primitive.py ${CMAKE_CURRENT_BINARY_DIR})

  add_test(smith python ${CMAKE_CURRENT_SOURCE_DIR}/smith.py ${CMAKE_CURRENT_BINARY_DIR})

  add_test(cellinvariants python ${CMAKE_CURRENT_SOURCE_DIR}/cellinvariants.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(space_group python ${CMAKE_CURRENT_SOURCE_DIR}/space_group.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(periodic_dnc python ${CMAKE_CURRENT_SOURCE_DIR}/periodic_dnc.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(neighbors python ${CMAKE_CURRENT_SOURCE_DIR}/neighbors.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(coordination_shells python ${CMAKE_CURRENT_SOURCE_DIR}/coordination_shells.py
           ${CMAKE_CURRENT_BINARY_DIR})
  add_test(splitconfigs python ${CMAKE_CURRENT_SOURCE_DIR}/splitconfigs.py
           ${CMAKE_CURRENT_BINARY_DIR})
  add_test(map_sites python ${CMAKE_CURRENT_SOURCE_DIR}/map_sites.py
           ${CMAKE_CURRENT_BINARY_DIR})

endif(NOT compile_debug_tests)
