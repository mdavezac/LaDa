set( LIBENUM_SOURCES 
     find_all_cells.cc label_exchange.cc numeric_type.cc
     remove_smith_equivalents.cc transform.cc translation.cc
     exceptions.h  find_all_cells.h  label_exchange.h  numeric_type.h
     remove_smith_equivalents.h  transform.h  translation.h )


add_library (enum STATIC ${LIBENUM_SOURCES})
target_link_libraries(enum crystal)
