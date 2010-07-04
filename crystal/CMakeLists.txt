set( LIBCRYSTAL_HEADERS 
     atom.h divide_and_conquer.h epi_structure.h fractional_cartesian.h layerdepth.h 
     read_structure.h structure.impl.h compare_sites.h divide_and_conquer.impl.h 
     fill_structure.h ideal_lattice.h neighbors.h smith.h symmetry_operator.h 
     confsplit.h enumerate.h fourier.h lattice.h read_poscar.h structure.h which_site.h )

set( LIBCRYSTAL_SOURCES 
     confsplit.cc fill_structure.cc lattice.cc
     read_structure.cc structure.cc epi_structure.cc ideal_lattice.cc 
     layerdepth.cc make_primitive.cc smith.cc symmetry_operator.cc )
    
add_library (crystal SHARED ${LIBCRYSTAL_SOURCES} ${LIBCRYSTAL_HEADERS})

target_link_libraries(crystal opt mathlib)

if(NOT noinstall_headers)
  install(FILES ${LIBCRYSTAL_HEADERS} DESTINATION include/LaDa/crystal)
endif(NOT noinstall_headers)

if(python_bindings)
  add_subdirectory(python)
endif(python_bindings)
