if(NOT compile_debug_tests)

  add_test(dftcrystal_rawgeom python ${CMAKE_CURRENT_SOURCE_DIR}/rawgeometry.py)
  add_test(dftcrystal_inputgeom python ${CMAKE_CURRENT_SOURCE_DIR}/inputgeom.py)

endif(NOT compile_debug_tests)
