if(NOT compile_debug_tests)

  add_test(dftcrystal_rawgeom python ${CMAKE_CURRENT_SOURCE_DIR}/rawgeometry.py)
  add_test(dftcrystal_rawbasis python ${CMAKE_CURRENT_SOURCE_DIR}/rawbasis.py)
  add_test(dftcrystal_keywords python ${CMAKE_CURRENT_SOURCE_DIR}/keywords.py)
  add_test(dftcrystal_optgeom python ${CMAKE_CURRENT_SOURCE_DIR}/optgeom.py)
  add_test(dftcrystal_hamiltonian python ${CMAKE_CURRENT_SOURCE_DIR}/hamiltonian.py)
  add_test(dftcrystal_electronic python ${CMAKE_CURRENT_SOURCE_DIR}/electronic.py)
  add_test(dftcrystal_geometry python ${CMAKE_CURRENT_SOURCE_DIR}/geometry.py)
  add_test(dftcrystal_inputgeom python ${CMAKE_CURRENT_SOURCE_DIR}/inputgeom.py)
  add_test(dftcrystal_pickle python ${CMAKE_CURRENT_SOURCE_DIR}/pickle.py)

endif(NOT compile_debug_tests)