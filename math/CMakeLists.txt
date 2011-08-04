set( LIBMATH_F90 rational_mathematics.f90 )
set( LIBMATH_HEADERS 
     compare_norms.h eigen.h fuzzy.h lambda.impl.h misc.h modifiers.h random.h serialize.h
     smith_normal_form.h traits.h symmetry_operator.h )
set( LIBMATH_SOURCES random.cc  smith_normal_form.cc )
    
add_library (mathlib SHARED ${LIBMATH_SOURCES} ${LIBMATH_F90} ${LIBMATH_HEADERS})
target_link_libraries(mathlib ${Boost_SERIALIZATION_LIBRARY} ${Boost_REGEX_LIBRARY})
set_target_properties(mathlib PROPERTIES OUTPUT_NAME lada_math)
install(TARGETS mathlib DESTINATION lib64)

include_directories(${tinyxml_INCLUDE_DIRS})
if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  if(MPI_LIBRARY) # crays make things complicated.
    target_link_libraries(mathlib  ${MPI_LIBRARY})
  endif(MPI_LIBRARY)
endif(MPI_FOUND)

if(NOT noinstall_headers)
  install(FILES ${LIBMATH_HEADERS} DESTINATION include/LaDa/math)
endif(NOT noinstall_headers)

if(python_bindings)
  add_subdirectory(python)
endif(python_bindings)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)

