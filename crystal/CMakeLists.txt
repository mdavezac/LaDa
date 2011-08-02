set( LIBCRYSTAL_HEADERS 
     atom.h structure.h structure_data.h add_atom.h is_container.h)
#    atom.h divide_and_conquer.h epi_structure.h fractional_cartesian.h layerdepth.h 
#    read_structure.h structure.impl.h compare_sites.h divide_and_conquer.impl.h 
#    fill_structure.h neighbors.h smith.h symmetry_operator.h 
#    confsplit.h enumerate.h fourier.h lattice.h read_poscar.h structure.h which_site.h 
#    periodic_dnc.h )

# set( LIBCRYSTAL_SOURCES  )
#    structure.cc )
#    confsplit.cc fill_structure.cc lattice.cc periodic_dnc.cc
#    read_structure.cc structure.cc epi_structure.cc neighbors.cc
#    layerdepth.cc make_primitive.cc smith.cc symmetry_operator.cc )
    
# add_library (crystal SHARED ${LIBCRYSTAL_HEADERS}) # ${LIBCRYSTAL_SOURCES} 
# set_target_properties(crystal PROPERTIES OUTPUT_NAME lada_crystal)
# target_link_libraries(crystal opt mathlib)
# install(TARGETS crystal DESTINATION lib64)

if(NOT noinstall_headers)
  install(FILES ${LIBCRYSTAL_HEADERS} DESTINATION include/LaDa/crystal)
endif(NOT noinstall_headers)

if(python_bindings)
  add_subdirectory(python)
  add_subdirectory(defects)
endif(python_bindings)
if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
