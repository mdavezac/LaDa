set( LIBMATH_F90 rational_mathematics.f90 )
set( LIBMATH_HEADERS 
     compare_norms.h eigen.h fuzzy.h lambda.impl.h misc.h modifiers.h random.h serialize.h
     smith_normal_form.h traits.h )
set( LIBMATH_SOURCES random.cc  smith_normal_form.cc )
    
add_library (mathlib STATIC ${LIBMATH_SOURCES} ${LIBMATH_F90} ${LIBMATH_HEADERS})
target_link_libraries(mathlib ${tinyxml_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_REGEX_LIBRARY})
set_target_properties(mathlib PROPERTIES OUTPUT_NAME math)
include_directories(${tinyxml_INCLUDE_DIRS})
if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  target_link_libraries(mathlib ${MPI_LIBRARY})
endif(MPI_FOUND)

if(NOT noinstall_headers)
  install(FILES ${LIBMATH_HEADERS} DESTINATION include/LaDa/math)
endif(NOT noinstall_headers)
