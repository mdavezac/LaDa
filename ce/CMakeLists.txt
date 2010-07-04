set( LIBCE_HEADERS 
     cluster.h constituent_strain.impl.h create_pairs.h find_pis.h fit.impl.h
     functional_builder.impl.h harmonic.impl.h  mlclusters.h regularization.h
     constituent_strain.h create_clusters.h drautz_diaz_ortiz.h  fit.h
     functional_builder.h  harmonic.h mlcluster.h polynome.h)
set( LIBCE_SOURCES 
     cluster.cc  create_clusters.cc  create_pairs.cc  fit.cc  harmonic.cc
     mlcluster.cc  mlclusters.cc  polynome.cc  regularization.cc  )



add_library (ce STATIC ${LIBCE_SOURCES} ${LIBCE_HEADERS})
target_link_libraries(ce crystal)
if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  target_link_libraries(ce ${MPI_LIBRARIES} ${Boost_MPI_LIBRARY})
endif (MPI_FOUND)

if(NOT disable_executables)
  set(CUBIC_CXXFLAGS "-D_CUBIC_CE_")
  set(TETRA_CXXFLAGS "-D_TETRAGONAL_CE_")
  set(ALL_LDFLAGS "")
  if(MPI_FOUND)
    set(CUBIC_CXXFLAGS ${CUBIC_CXXFLAGS} ${MPI_COMPILE_FLAGS})
    set(TETRA_CXXFLAGS ${CUBIC_CXXFLAGS} ${MPI_COMPILE_FLAGS})
    set(ALL_LDFLAGS ${ALL_LDFLAGS} ${MPI_LINK_FLAGS})
  endif(MPI_FOUND)

  add_executable(cubic main.cc)
  set_target_properties(cubic PROPERTIES LINK_FLAGS "${ALL_LDFLAGS}")
  set_target_properties(cubic PROPERTIES COMPILE_FLAGS "${CUBIC_CXXFLAGS}")
  target_link_libraries(cubic ce ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY})
  
  add_executable(tetra main.cc)
  set_target_properties(tetra PROPERTIES LINK_FLAGS "${ALL_LDFLAGS}")
  set_target_properties(tetra PROPERTIES COMPILE_FLAGS "${TETRA_CXXFLAGS}")
  target_link_libraries(tetra ce ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY})

  add_executable(regularized_cubic rmain.cc)
  set_target_properties(regularized_cubic PROPERTIES LINK_FLAGS "${ALL_LDFLAGS}")
  set_target_properties(regularized_cubic PROPERTIES COMPILE_FLAGS "${CUBIC_CXXFLAGS}")
  target_link_libraries(regularized_cubic ce minimizer ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY})
  add_executable(regularized_tetra rmain.cc)
  set_target_properties(regularized_tetra PROPERTIES LINK_FLAGS "${ALL_LDFLAGS}")
  set_target_properties(regularized_tetra PROPERTIES COMPILE_FLAGS "${TETRA_CXXFLAGS}")
  target_link_libraries(regularized_tetra ce minimizer ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY})

  install(TARGETS cubic DESTINATION bin)
  install(TARGETS tetra DESTINATION bin)
  install(TARGETS regularized_cubic DESTINATION bin)
  install(TARGETS regularized_tetra DESTINATION bin)
endif(NOT disable_executables)

if(NOT noinstall_headers)
  install(FILES ${LIBCE_HEADERS} DESTINATION include/LaDa/ce)
endif(NOT noinstall_headers)

