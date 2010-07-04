set( LIBENUM_HEADERS 
     exceptions.h  find_all_cells.h  label_exchange.h  numeric_type.h
     remove_smith_equivalents.h  transform.h  translation.h )
set( LIBENUM_SOURCES 
     find_all_cells.cc label_exchange.cc numeric_type.cc
     remove_smith_equivalents.cc transform.cc translation.cc )


add_library (enum SHARED ${LIBENUM_SOURCES} ${LIBENUM_HEADERS})
target_link_libraries(enum crystal)
if(NOT noinstall_headers)
  install(FILES ${LIBENUM_HEADERS} DESTINATION include/LaDa/enum)
endif(NOT noinstall_headers)
