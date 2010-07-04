set( LIBOPT_F90  fortran_redirect.f90)
set( LIBOPT_SOURCES  convex_hull.cc  errors.cc leave_many_out.cc path.cc tinyxml.cc )
set( LIBOPT_HEADERS 
     algorithms.h convex_hull.h debug.h fitness_function.h function_base.h
     indirection.h monome.h path.h pow.h redirect.h tuple_io.h types.h
     bpo_macros.h convex_hull.impl.h errors.h function_functors.h
     leave_many_out.h ndim_iterator.h polynome.h propagate_std_vector.h
     tinyxml.h tuple_serialize.h va_function.h
   )
    
add_library (opt STATIC ${LIBOPT_HEADERS} ${LIBOPT_SOURCES} ${LIBOPT_F90})
target_link_libraries(opt ${tinyxml_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_REGEX_LIBRARY}
                          ${Boost_SERIALIZATION_LIBRARY})
include_directories(${tinyxml_INCLUDE_DIRS})
if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  target_link_libraries(opt ${MPI_LIBRARY})
endif(MPI_FOUND)

if(NOT noinstall_headers)
  install(FILES ${LIBOPT_HEADERS} DESTINATION include/LaDa/opt)
endif(NOT noinstall_headers)

if(python_bindings)
  add_subdirectory(python)
endif(python_bindings)
