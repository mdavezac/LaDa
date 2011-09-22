if(NOT compile_debug_tests)
  add_test(pyset python ${CMAKE_CURRENT_SOURCE_DIR}/pyset.py)
  add_test(atoms python ${CMAKE_CURRENT_BINARY_DIR}/atoms.py)
  add_test(pickle_atom python ${CMAKE_CURRENT_BINARY_DIR}/pickle_atom.py)
  add_test(structure python ${CMAKE_CURRENT_BINARY_DIR}/structure.py)
  add_test(structure_contains python ${CMAKE_CURRENT_BINARY_DIR}/structure_contains.py)
endif(NOT compile_debug_tests)
  add_test(pickle_structure python ${CMAKE_CURRENT_BINARY_DIR}/pickle_structure.py)

  
