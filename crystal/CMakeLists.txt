set(LIBCRYSTAL_PYSOURCES __init__.py read.py write.py ABX.py A2BX4.py bravais.py binary.py)

include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})

add_library (pycrystal SHARED module.cc)
target_link_libraries(pycrystal ${PYTHON_LIBRARY} mathlib)
set_target_properties(pycrystal PROPERTIES OUTPUT_NAME cppwrappers)
set_target_properties(pycrystal PROPERTIES PREFIX "")
install(TARGETS pycrystal DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal)
install(FILES ${LIBCRYSTAL_PYSOURCES} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
