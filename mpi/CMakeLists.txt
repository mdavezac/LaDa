set(LIBMPI_HEADERS mpi_object.h macros.h)
if(NOT noinstall_headers)
  install(FILES ${LIBMPI_HEADERS} DESTINATION include/LaDa/mpi)
endif(NOT noinstall_headers)
