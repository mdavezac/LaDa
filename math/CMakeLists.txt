set(PYMATH_HEADERS fuzzy.h gruber.h eigen.h smith_normal_form.h exceptions.h
    misc.h)
set(PYMATH_SOURCES methods.hpp module.cc gruber.cc noopt.cc smith_normal_form.cc)

include_directories(${PYTHON_INCLUDE_DIR} ${NUMPY_INCLUDES})
add_library (pymath SHARED ${PYMATH_HEADERS} ${PYMATH_SOURCES})
SET_SOURCE_FILES_PROPERTIES(noopt.cc PROPERTIES COMPILE_FLAGS -O0)
set_target_properties(pymath PROPERTIES OUTPUT_NAME math)
set_target_properties(pymath PROPERTIES PREFIX "")
target_link_libraries(pymath ${PYTHON_LIBRARY})
install(TARGETS pymath DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)

