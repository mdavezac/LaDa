set(LIBCRYSTAL_HEADERS atom.h atom_base.h structure_data.h)
set(LIBCRYSTAL_PYHEADERS atom_cdi.hpp atom_getset.hpp structure_cdi.hpp structure_getset.hpp)
set(LIBCRYSTAL_SOURCES atom_base.cc structure_data.cc )

include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
add_library (crystal SHARED ${LIBCRYSTAL_HEADERS} ${LIBCRYSTAL_PYHEADERS} ${LIBCRYSTAL_SOURCES})
target_link_libraries(crystal ${PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARY} mathlib)
install(TARGETS crystal DESTINATION lib64)

add_library (pycrystal SHARED module.cc)
target_link_libraries(pycrystal ${PYTHON_LIBRARY} crystal)
set_target_properties(pycrystal PROPERTIES OUTPUT_NAME cppwrappers)
set_target_properties(pycrystal PROPERTIES PREFIX "")
install(TARGETS pycrystal DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal)
install(FILES python/__init__.py DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal)

if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
endif(MPI_FOUND)
if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)

# set( LIBCRYSTAL_HEADERS 
#      atom.h atom_base.h atom_freeze.h structure.h structure_data.h add_atom.h is_container.h
#      utilities.h compare_sites.h space_group.h neighbors.h map_sites.h 
#      periodic_dnc.h cast.h)
# #    atom.h divide_and_conquer.h epi_structure.h fractional_cartesian.h layerdepth.h 
# #    read_structure.h structure.impl.h compare_sites.h divide_and_conquer.impl.h 
# #    fill_structure.h neighbors.h smith.h symmetry_operator.h 
# #    confsplit.h enumerate.h fourier.h lattice.h read_poscar.h structure.h which_site.h 
# #    periodic_dnc.h )

# set( LIBCRYSTAL_SOURCES utilities.cc space_group.cc periodic_dnc.cc cast.cc)
# #    structure.cc )
# #    confsplit.cc fill_structure.cc lattice.cc periodic_dnc.cc
# #    read_structure.cc structure.cc epi_structure.cc neighbors.cc
# #    layerdepth.cc make_primitive.cc smith.cc symmetry_operator.cc )
#     
# add_library (crystal SHARED ${LIBCRYSTAL_HEADERS} ${LIBCRYSTAL_SOURCES})
# set_target_properties(crystal PROPERTIES OUTPUT_NAME lada_crystal)
# if(loadnsave)
#   target_link_libraries(crystal lns)
# endif(loadnsave)
# if(python_bindings)
#   include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
#   target_link_libraries(crystal ${Boost_PYTHON_LIBRARIES} ${PYTHON_LIBRARY})
# endif(python_bindings)
# target_link_libraries(crystal mathlib)
# install(TARGETS crystal DESTINATION lib64)

# if(NOT noinstall_headers)
#   install(FILES ${LIBCRYSTAL_HEADERS} DESTINATION include/LaDa/crystal)
# endif(NOT noinstall_headers)

# if(python_bindings)
#   include_directories(${PYTHON_INCLUDE_DIRS})
#   add_subdirectory(python)
#   add_subdirectory(defects)
# endif(python_bindings)
# if(compile_tests OR compile_debug_tests)
#   add_subdirectory(tests)
# endif(compile_tests OR compile_debug_tests)
