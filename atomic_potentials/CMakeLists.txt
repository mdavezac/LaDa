set( _LIBCOLLAPSE_HEADERS 
     collapse.h fitting_set.h values.h values.rank_iterator.h values.rep_iterator.h 
     values.str_iterator.h )
set( _LIBCOLLAPSE_SOURCES  collapse.cc fitting_set.cc values.cc )
set( LIBATOMPOT_HEADERS 
     bases.h bases.impl.h coord_range.rank_range.h  functions.h
     functions.iterator.h numeric_types.h representation.h  separable.h
     sum_of_separables.coord_range.h  sum_of_separables.h )
set( LIBATOMPOT_SOURCES functions.cc  representation.cc )
foreach(ATOMPOTFILE ${_LIBCOLLAPSE_SOURCES})
  set(LIBCOLLAPS_SOURCES ${LIBCOLLAPSE_SOURCES} collapse/${ATOMPOTFILE})
endforeach(ATOMPOTFILE ${_LIBCOLLAPSE_SOURCES})
foreach(ATOMPOTFILE ${_LIBCOLLAPSE_HEADERS})
  set(LIBCOLLAPSE_HEADERS ${LIBCOLLAPSE_HEADERS} collapse/${ATOMPOTFILE})
endforeach(ATOMPOTFILE ${_LIBCOLLAPSE_HEADERS})

add_library (atompot STATIC ${LIBATOMPOT_SOURCES} ${LIBATOM_HEADERS} ${LIBCOLLAPSE_HEADERS} ${LIBCOLLAPSE_SOURCES})
target_link_libraries(atompot crystal)

if(NOT noinstall_headers)
  install(FILES ${LIBATOMPOT_HEADERS} DESTINATION include/LaDa/atomic_potentials)
  install(FILES ${LIBCOLLAPSE_HEADERS} DESTINATION include/LaDa/atomic_potentials/collapse)
endif(NOT noinstall_headers)

