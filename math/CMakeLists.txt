set(LIBMATH_HEADERS fuzzy.h gruber.h eigen.h smith_normal_form.h exceptions.h
                    misc.h )
set(LIBMATH_SOURCES gruber.cc noopt.cc smith_normal_form.cc )
set(PYMATH_SOURCES methods.hpp module.cc gruber.cc noopt.cc smith_normal_form.cc)

add_library (mathlib SHARED ${LIBMATH_HEADERS} ${LIBMATH_SOURCES})
SET_SOURCE_FILES_PROPERTIES(noopt.cc PROPERTIES COMPILE_FLAGS -O0)
target_link_libraries(mathlib ${PYTHON_LIBRARY})
install(TARGETS mathlib DESTINATION lib64)

include_directories(${PYTHON_INCLUDE_DIR} ${NUMPY_INCLUDES})
add_library (pymath SHARED ${PYMATH_SOURCES})
set_target_properties(pymath PROPERTIES OUTPUT_NAME math)
set_target_properties(pymath PROPERTIES PREFIX "")
target_link_libraries(pymath mathlib)
install(TARGETS pymath DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)

