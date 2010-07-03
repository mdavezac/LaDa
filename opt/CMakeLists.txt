
set( LIBOPT_SOURCES 
     convex_hull.cc  errors.cc  fortran_dirio.cc  leave_many_out.cc  path.cc tinyxml.cc 
     algorithms.h convex_hull.h debug.h fitness_function.h function_base.h indirection.h 
     monome.h path.h pow.h redirect.h tuple_io.h types.h bpo_macros.h convex_hull.impl.h 
     errors.h fortran_dirio.h function_functors.h leave_many_out.h ndim_iterator.h polynome.h
     propagate_std_vector.h tinyxml.h tuple_serialize.h va_function.h
   )
    
add_library (opt SHARED ${LIBOPT_SOURCES})
