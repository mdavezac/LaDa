set( LIBCOLLAPSE_SOURCES 
     collapse.cc fitting_set.cc values.cc
     collapse.h fitting_set.h values.h values.rank_iterator.h values.rep_iterator.h 
     values.str_iterator.h )
set( LIBATOMPOT_SOURCES 
     functions.cc  representation.cc
     bases.h bases.impl.h coord_range.rank_range.h  functions.h
     functions.iterator.h numeric_types.h representation.h  separable.h
     sum_of_separables.coord_range.h  sum_of_separables.h )
foreach(ATOMPOTFILE ${LIBCOLLAPSE_SOURCES})
  set(LIBATOMPOT_SOURCES ${LIBATOMPOT_SOURCES} collapse/${ATOMPOTFILE})
endforeach(ATOMPOTFILE ${LIBCOLLAPSE_SOURCES})



add_library (atompot STATIC ${LIBATOMPOT_SOURCES})
target_link_libraries(atompot crystal)

