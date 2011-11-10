if(NOT compile_debug_tests)
  add_test(atoms python ${CMAKE_CURRENT_BINARY_DIR}/atoms.py)
  add_test(pickle_atom python ${CMAKE_CURRENT_BINARY_DIR}/pickle_atom.py)
  add_test(structure python ${CMAKE_CURRENT_BINARY_DIR}/structure.py)
  add_test(structure_contains python ${CMAKE_CURRENT_BINARY_DIR}/structure_contains.py)
  add_test(pickle_structure python ${CMAKE_CURRENT_BINARY_DIR}/pickle_structure.py)
endif(NOT compile_debug_tests)
  add_test(pyset python ${CMAKE_CURRENT_SOURCE_DIR}/pyset.py)

set(TYPES AtomStr AtomSet)

# foreach(TYPE ${TYPES})
#   if(${TYPE} STREQUAL "AtomStr")
#     set(LADA_CTYPE "std::string")
#   elseif(${TYPE} STREQUAL "AtomVec")
#     set(LADA_CTYPE "std::vector<std::string>")
#   elseif(${TYPE} STREQUAL "AtomSet")
#     set(LADA_CTYPE "std::set<std::string>")
#   endif(${TYPE} STREQUAL "AtomStr")

#   python_add_module (atom_self_${TYPE} "" "atom_self.cc")
#   set_property( TARGET atom_self_${TYPE} PROPERTY COMPILE_DEFINITIONS
#                 LADA_TYPE=${LADA_CTYPE} LADA_MODULE=atom_self_${TYPE})
#   mark_as_advanced(PYTHON_ENABLE_MODULE_atom_self_${TYPE})
#   mark_as_advanced(PYTHON_MODULE_atom_self_${TYPE}_BUILD_SHARED)
#   set_target_properties(atom_self_${TYPE} PROPERTIES OUTPUT_NAME atom_self_${TYPE})
#   target_link_libraries(atom_self_${TYPE} ${Boost_PYTHON_LIBRARY})
#   configure_file (
#    "${CMAKE_CURRENT_SOURCE_DIR}/atom_self_py.py"
#    "${CMAKE_CURRENT_BINARY_DIR}/atom_self_${TYPE}.py"
#   )
#   add_test(atom_self_${TYPE} python ${CMAKE_CURRENT_BINARY_DIR}/atom_self_${TYPE}.py)
# endforeach(TYPE ${TYPES})
