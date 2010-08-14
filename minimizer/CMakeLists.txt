set( LIBMINIMIZER_F90 df1dim.f90  frprmn_fortran.f90  linmin.f90 )
set( LIBMINIMIZER_HEADERS 
     allsq.h  cgs.h  decoupled.h  frprmn.h  function_wrapper.h
     interpolated_gradient.h variant.h )
set( LIBMINIMIZER_SOURCES cgs.cc decoupled.cc frprmn.cc )

if (GSL_FOUND)
  include_directories(${gsl_INCLUDE_DIRS})
  set( LIBMINIMIZER_HEADERS  ${LIBMINIMIZER_HEADERS} 
       gsl.h  gsl_lsq.h  gsl_mins.h  gsl_nllsq.h  gsl_simplex.h)
  set( LIBMINIMIZER_SOURCES  ${LIBMINIMIZER_SOURCES} 
       gsl.cc gsl_lsq.cc gsl_mins.cc gsl_nllsq.cc gsl_simplex.cc )
endif (GSL_FOUND)

if (MINUIT2_FOUND)
  include_directories(${minuit2_INCLUDE_DIRS})
  set(LIBMINIMIZER_SOURCES  ${LIBMINIMIZER_SOURCES} minuit2.cc )
  set(LIBMINIMIZER_HEADERS  ${LIBMINIMIZER_HEADERS} minuit2.h )
endif (MINUIT2_FOUND)

    
add_library (minimizer SHARED ${LIBMINIMIZER_SOURCES} ${LIBMINIMIZER_HEADERS} ${LIBMINIMIZER_F90})
set_target_properties(minimizer PROPERTIES OUTPUT_NAME lada_minimizer)
target_link_libraries(minimizer opt)
install(TARGETS minimizer DESTINATION lib64)
if(NOT noinstall_headers)
  install(FILES ${LIBMINIMIZER_HEADERS} DESTINATION include/LaDa/minimizer)
endif(NOT noinstall_headers)

if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  if(MPI_LIBRARY) # crays make things complicated.
    target_link_libraries(minimizer  ${MPI_LIBRARY})
  endif(MPI_LIBRARY)
endif (MPI_FOUND)
if (GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  target_link_libraries(minimizer ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})
endif(GSL_FOUND)
if (MINUIT2_FOUND)
  include_directories(${minuit2_INCLUDE_DIRS})
  target_link_libraries(minimizer ${minuit2_LIBRARY})
endif(MINUIT2_FOUND)

if(python_bindings)
  add_subdirectory(python)
endif(python_bindings)
