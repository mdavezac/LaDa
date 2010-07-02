
set( LIBCRYSTAL_SOURCES \
     confsplit.cc fill_structure.cc lattice.cc main.cc \
     read_structure.cc structure.cc epi_structure.cc ideal_lattice.cc \
     layerdepth.cc make_primitive.cc smith.cc symmetry_operator.cc  )
add_library (crystal SHARED ${LIBCRYSTAL_SOURCES})
