set( LIBMINIMIZER_SOURCES 
     allsq.h  cgs.h  decoupled.h  frprmn.h  function_wrapper.h interpolated_gradient.h variant.h
     cgs.cc decoupled.cc frprmn.cc 
     df1dim.f90  frprmn_fortran.f90  linmin.f90
   )

if (GSL_FOUND)
  set( LIBMINIMIZER_SOURCES  ${LIBMINIMIZER_SOURCES} 
       gsl.h  gsl_lsq.h  gsl_mins.h  gsl_nllsq.h  gsl_simplex.h
       gsl.cc gsl_lsq.cc gsl_mins.cc gsl_nllsq.cc gsl_simplex.cc 
     )
endif (GSL_FOUND)

if (MINUIT2_FOUND)
  set(LIBMINIMIZER_SOURCES  ${LIBMINIMIZER_SOURCES} minuit2.h minuit2.cc)
endif (MINUIT2_FOUND)

    
add_library (minimizer SHARED ${LIBMINIMIZER_SOURCES})
target_link_libraries(minimizer opt)
target_link_libraries(minimizer ${Boost_python_LIBRARY})

if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  target_link_libraries(minimizer ${MPI_LIBRARY})
endif (MPI_FOUND)
if (GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  target_link_libraries(minimizer ${GSL_LIBRARY})
endif(GSL_FOUND)
if (MINUIT2_FOUND)
  include_directories(${MINUIT2_INCLUDE_DIRS})
  target_link_libraries(minimizer ${MINUIT2_LIBRARY})
endif(MINUIT2_FOUND)
