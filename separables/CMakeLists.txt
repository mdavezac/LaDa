set( LIBFXSEPS_HEADERS 
     ce_as_collapse.h collapse.h colmappings.h colpolicy.impl.h
     functional.impl.h many.impl.h methods.impl.h mixedfunctional.impl.h
     mixed.impl.h sepmappings.h ce_as_collapse.impl.h collapse.impl.h
     colpolicy.h functional.h many.h methods.h mixedfunctional.h mixed.h
     prepare.h )
set( LIBFXSEPS_SOURCES prepare.cc )


add_library (fxseps SHARED ${LIBFXSEPS_SOURCES} ${LIBFXSEPS_HEADERS})
target_link_libraries(fxseps ce)
if(NOT noinstall_headers)
  install(FILES ${LIBFXSEPS_HEADERS} DESTINATION include/LaDa/separables)
endif(NOT noinstall_headers)

if(NOT disable_executables)
  set(CUBIC_CXXFLAGS "-D_CUBIC_CE_")
  set(TETRA_CXXFLAGS "-D_TETRAGONAL_CE_")
  set(ALL_LDFLAGS "")
  if(MPI_FOUND)
    set(CUBIC_CXXFLAGS ${CUBIC_CXXFLAGS} ${MPI_COMPILE_FLAGS})
    set(TETRA_CXXFLAGS ${CUBIC_CXXFLAGS} ${MPI_COMPILE_FLAGS})
    set(ALL_LDFLAGS ${ALL_LDFLAGS} ${MPI_LINK_FLAGS})
    include_directories(${MPI_INCLUDE_PATH})
  endif(MPI_FOUND)


  add_executable(seps_cubic main.cc)
  set_target_properties(seps_cubic PROPERTIES LINK_FLAGS "${ALL_LDFLAGS}")
  set_target_properties(seps_cubic PROPERTIES COMPILE_FLAGS "${CUBIC_CXXFLAGS}")
  target_link_libraries(seps_cubic fxseps minimizer ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY})

  add_executable(seps_tetra main.cc)
  set_target_properties(seps_tetra PROPERTIES LINK_FLAGS "${ALL_LDFLAGS}")
  set_target_properties(seps_tetra PROPERTIES COMPILE_FLAGS "${TETRA_CXXFLAGS}")
  target_link_libraries(seps_tetra fxseps minimizer ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY})

  install(TARGETS seps_cubic DESTINATION bin)
  install(TARGETS seps_tetra DESTINATION bin)
endif(NOT disable_executables)

