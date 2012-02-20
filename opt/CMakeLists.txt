set( LIBOPT_F90  fortran_redirect.F90)
set( LIBOPT_SOURCES  convex_hull.cc  path.cc )
# errors.cc
set( LIBOPT_HEADERS 
     algorithms.h convex_hull.h debug.h fitness_function.h function_base.h
     monome.h path.h pow.h redirect.h tuple_io.h types.h
     convex_hull.impl.h function_functors.h
     ndim_iterator.h polynome.h propagate_std_vector.h
     tuple_serialize.h va_function.h
   ) 
# errors.h bpo_macros.h leave_many_out.h leave_many_out.cc tinyxml.h indirection.h 
    
add_library (opt SHARED ${LIBOPT_HEADERS} ${LIBOPT_SOURCES} ${LIBOPT_F90})
set_target_properties(opt PROPERTIES OUTPUT_NAME lada_opt)

target_link_libraries(opt mathlib ${tinyxml_LIBRARY}
                          ${Boost_SYSTEM_LIBRARY} 
                          ${Boost_FILESYSTEM_LIBRARY} 
                          ${Boost_REGEX_LIBRARY}
                          ${Boost_SERIALIZATION_LIBRARY})
include_directories(${tinyxml_INCLUDE_DIRS})
if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  target_link_libraries(opt ${Boost_MPI_LIBRARY})
  if(MPI_LIBRARY) # crays make things complicated.
    target_link_libraries(opt  ${MPI_LIBRARY})
  endif(MPI_LIBRARY)
endif(MPI_FOUND)

if(NOT noinstall_headers)
  install(FILES ${LIBOPT_HEADERS} DESTINATION include/LaDa/opt)
endif(NOT noinstall_headers)

install(TARGETS opt DESTINATION lib64)

#  if(python_bindings)
#    add_subdirectory(python)
#  endif(python_bindings)
