set(PYMATH_HEADERS methods.hpp fuzzy.h)
set(PYMATH_SOURCES module.cc)

include_directories(${PYTHON_INCLUDE_DIR} ${NUMPY_INCLUDES})
add_library (pymath SHARED ${PYMATH_HEADERS} ${PYMATH_SOURCES})
set_target_properties(pymath PROPERTIES OUTPUT_NAME math)
set_target_properties(pymath PROPERTIES PREFIX "")
install(TARGETS pymath DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada)
# set( LIBMATH_HEADERS 
#      compare_norms.h eigen.h fuzzy.h lambda.impl.h misc.h modifiers.h serialize.h
#      smith_normal_form.h symmetry_operator.h gruber.h quantity.h )
# set( LIBMATH_SOURCES smith_normal_form.cc gruber.cc noopt.cc quantity.cc)
#     
# add_library (mathlib SHARED ${LIBMATH_SOURCES} ${LIBMATH_HEADERS})
# target_link_libraries(mathlib ${PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARY} 
#                               ${Boost_SERIALIZATION_LIBRARY} ${Boost_REGEX_LIBRARY})
# set_target_properties(mathlib PROPERTIES OUTPUT_NAME lada_math)
# install(TARGETS mathlib DESTINATION lib64)
# SET_SOURCE_FILES_PROPERTIES(noopt.cc PROPERTIES COMPILE_FLAGS -O0)

# include_directories(${PYTHON_INCLUDE_DIRS})
# if(MPI_FOUND)
#   include_directories(${MPI_INCLUDE_PATH})
#   if(MPI_LIBRARY) # crays make things complicated.
#     target_link_libraries(mathlib  ${MPI_LIBRARY})
#   endif(MPI_LIBRARY)
# endif(MPI_FOUND)

# if(NOT noinstall_headers)
#   install(FILES ${LIBMATH_HEADERS} DESTINATION include/LaDa/math)
# endif(NOT noinstall_headers)

# if(compile_tests OR compile_debug_tests)
#   add_subdirectory(tests)
# endif(compile_tests OR compile_debug_tests)

# add_subdirectory(python)
