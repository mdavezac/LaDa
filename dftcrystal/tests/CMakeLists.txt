if(NOT compile_debug_tests)

  add_test(dftcrystal_rawgeom python ${CMAKE_CURRENT_SOURCE_DIR}/rawgeometry.py ${CMAKE_CURRENT_SOURCE_DIR})

endif(NOT compile_debug_tests)
