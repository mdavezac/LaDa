set( LIBVFF_HEADERS 
     atomic_center.h atomic_functional.h functional.h layered.h va.h va.impl.h vff.h )
set( LIBVFF_SOURCES 
     atomic_center.cc atomic_functional.cc build_tree_smith.cc build_tree_sort.cc
     build_tree_sort_dnc.cc functional.cc initialize_centers.cc layered.cc vff.cc )


add_library (vff STATIC ${LIBVFF_SOURCES} ${LIBVFF_HEADERS})
target_link_libraries(vff crystal minimizer)
if(NOT noinstall_headers)
  install(FILES ${LIBVFF_HEADERS} DESTINATION include/LaDa/vff)
endif(NOT noinstall_headers)

if(NOT disable_executables)
  add_executable(vffexec main.cc)
  set_target_properties(vffexec PROPERTIES OUTPUT_NAME vff)
  if(MPI_FOUND)
    set_target_properties(vffexec PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    set_target_properties(vffexec PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif(MPI_FOUND)
  target_link_libraries(vffexec vff ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY})
  install(TARGETS vffexec DESTINATION bin)
endif(NOT disable_executables)


if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  target_link_libraries(vff ${MPI_LIBRARIES} ${Boost_MPI_LIBRARY})
endif (MPI_FOUND)
if (GSL_FOUND)
  include_directories(${gsl_INCLUDE_DIRS})
endif(GSL_FOUND)
if (MINUIT2_FOUND)
  include_directories(${minuit2_INCLUDE_DIRS})
endif(MINUIT2_FOUND)

