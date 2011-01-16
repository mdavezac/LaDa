set( LIBVFF_HEADERS 
     data.h exceptions.h atomic_center.h functional.h layered.h va.h va.impl.h vff.h )
set( LIBVFF_SOURCES 
     atomic_center.cc build_tree_smith.cc build_tree_sort.cc
     build_tree_sort_dnc.cc functional.cc initialize_centers.cc layered.cc vff.cc )


add_library (vff SHARED ${LIBVFF_SOURCES} ${LIBVFF_HEADERS})
set_target_properties(vff PROPERTIES OUTPUT_NAME lada_vff)
target_link_libraries(vff crystal minimizer)
install(TARGETS vff DESTINATION lib64)

if(NOT noinstall_headers)
  install(FILES ${LIBVFF_HEADERS} DESTINATION include/LaDa/vff)
endif(NOT noinstall_headers)

if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  target_link_libraries(vff ${MPI_LIBRARY} ${Boost_MPI_LIBRARY})
endif (MPI_FOUND)
if (GSL_FOUND)
  include_directories(${gsl_INCLUDE_DIRS})
endif(GSL_FOUND)
if (MINUIT2_FOUND)
  include_directories(${minuit2_INCLUDE_DIRS})
endif(MINUIT2_FOUND)

if(python_bindings)
  add_subdirectory(python)
endif(python_bindings)

