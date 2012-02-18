set(LIBCRYSTAL_HEADERS atom.h atom_base.h structure.h structure_data.h smith.h smith_data.h
                       supercell.h primitive.h space_group.h equivalent_structures.h
                       periodic_dnc.h neighbors.h coordination_shells.h confsplit.h )
set(LIBCRYSTAL_PYHEADERS atom_cdi.hpp atom_getset.hpp
                         structure_cdi.hpp structure_getset.hpp structure_sequence.hpp  
                         smith_cdi.hpp smith_get.hpp smith_members.hpp smith_macro.hpp)
set(LIBCRYSTAL_SOURCES atom_base.cc structure_data.cc structure.cc utilities.cc primitive.cc 
                       space_group.cc equivalent_structures.cc smith_data.cc periodic_dnc.cc
                       neighbors.cc coordination_shells.cc confsplit.cc )

include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
add_library (crystal SHARED ${LIBCRYSTAL_HEADERS} ${LIBCRYSTAL_PYHEADERS} ${LIBCRYSTAL_SOURCES})
target_link_libraries(crystal ${PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARY} mathlib)
install(TARGETS crystal DESTINATION lib64)

add_library (pycrystal SHARED module.cc)
target_link_libraries(pycrystal ${PYTHON_LIBRARY} crystal)
set_target_properties(pycrystal PROPERTIES OUTPUT_NAME cppwrappers)
set_target_properties(pycrystal PROPERTIES PREFIX "")
install(TARGETS pycrystal DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal)
install(FILES __init__.py DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
