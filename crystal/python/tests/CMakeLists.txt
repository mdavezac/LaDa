if(NOT compile_debug_tests)
# add_test(structure python ${CMAKE_CURRENT_BINARY_DIR}/structure.py)
# add_test(structure_contains python ${CMAKE_CURRENT_BINARY_DIR}/structure_contains.py)
# add_test(pickle_structure python ${CMAKE_CURRENT_BINARY_DIR}/pickle_structure.py)

  set(TYPES 0 1)
  foreach(TYPE ${TYPES})
    if(${TYPE} EQUAL 0)
      set(PYTYPE "AtomStr")
    elseif(${TYPE} EQUAL 1)
      set(PYTYPE "AtomSequence")
    endif(${TYPE} EQUAL 0)
    python_add_module (atom_self_${TYPE} "" "atom_self.cc")
    set_property(TARGET atom_self_${TYPE} PROPERTY COMPILE_DEFINITIONS LADA_TYPE=${TYPE})
    mark_as_advanced(PYTHON_ENABLE_MODULE_atom_self_${TYPE})
    mark_as_advanced(PYTHON_MODULE_atom_self_${TYPE}_BUILD_SHARED)
    set_target_properties(atom_self_${TYPE} PROPERTIES OUTPUT_NAME atom_self_${TYPE})
    target_link_libraries(atom_self_${TYPE} ${Boost_PYTHON_LIBRARY} pyatom)
    configure_file (
     "${CMAKE_CURRENT_SOURCE_DIR}/atom_self_py.py"
     "${CMAKE_CURRENT_BINARY_DIR}/atom_self_${TYPE}.py"
    )
    add_test(atom_self_${TYPE} python ${CMAKE_CURRENT_BINARY_DIR}/atom_self_${TYPE}.py)

    configure_file (
     "${CMAKE_CURRENT_SOURCE_DIR}/atom_memory.in.py"
     "${CMAKE_CURRENT_BINARY_DIR}/atom_memory_${TYPE}.py"
    )
    add_test(memory_${TYPE} python ${CMAKE_CURRENT_BINARY_DIR}/memstr_${TYPE}.py)
  endforeach(TYPE ${TYPES})
  add_test(pysequence python ${CMAKE_CURRENT_SOURCE_DIR}/pysequence.py)
  add_test(atomstr python ${CMAKE_CURRENT_SOURCE_DIR}/atomstr.py)
  add_test(atomseq python ${CMAKE_CURRENT_SOURCE_DIR}/atomseq.py)

endif(NOT compile_debug_tests)
