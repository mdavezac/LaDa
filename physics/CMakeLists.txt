set(LIBPHYSICS_HEADERS physics.h)
if(NOT noinstall_headers)
  install(FILES ${LIBPHYSICS_HEADERS} DESTINATION include/LaDa/atomic_potentials/physics)
endif(NOT noinstall_headers)
