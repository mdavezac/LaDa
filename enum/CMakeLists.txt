set (PYCE_PY __init__.py)

install(FILES ${PYCE_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/enum)

# set( LIBENUM_HEADERS 
#      exceptions.h  find_all_cells.h  label_exchange.h  numeric_type.h
#      remove_smith_equivalents.h  transform.h  translation.h )
# set( LIBENUM_SOURCES 
#      find_all_cells.cc label_exchange.cc numeric_type.cc
#      remove_smith_equivalents.cc transform.cc translation.cc )


# add_library (enum SHARED ${LIBENUM_SOURCES} ${LIBENUM_HEADERS})
# set_target_properties(enum PROPERTIES OUTPUT_NAME lada_enum)
# target_link_libraries(enum crystal)
# install(TARGETS enum DESTINATION lib64)

# if(NOT noinstall_headers)
#   install(FILES ${LIBENUM_HEADERS} DESTINATION include/LaDa/enum)
# endif(NOT noinstall_headers)

# if(python_bindings)
#   add_subdirectory(python)
# endif(python_bindings)

