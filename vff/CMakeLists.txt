set( LIBVFF_SOURCES 
     atomic_center.h atomic_functional.h functional.h layered.h va.h va.impl.h vff.h
     atomic_center.cc atomic_functional.cc build_tree_smith.cc build_tree_sort.cc
     build_tree_sort_dnc.cc functional.cc initialize_centers.cc layered.cc vff.cc )


add_library (vff SHARED ${LIBVFF_SOURCES})
target_link_libraries(vff crystal minimizer)
target_link_libraries(vff ${Boost_python_LIBRARY})

if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  target_link_libraries(vff ${MPI_LIBRARY})
endif (MPI_FOUND)
# if (GSL_FOUND)
#   include_directories(${GSL_INCLUDE_DIRS})
#   target_link_libraries(vff ${GSL_LIBRARY})
# endif(GSL_FOUND)
# if (MINUIT2_FOUND)
#   include_directories(${MINUIT2_INCLUDE_DIRS})
#   target_link_libraries(vff ${MINUIT2_LIBRARY})
# endif(MINUIT2_FOUND)
