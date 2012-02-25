if(NOT compile_debug_tests)
  
  add_test(algo python ${CMAKE_CURRENT_SOURCE_DIR}/algo.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(choices python ${CMAKE_CURRENT_SOURCE_DIR}/choices.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(magmom python ${CMAKE_CURRENT_SOURCE_DIR}/magmom.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(nelect python ${CMAKE_CURRENT_SOURCE_DIR}/nelect.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(npar python ${CMAKE_CURRENT_SOURCE_DIR}/npar.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(others python ${CMAKE_CURRENT_SOURCE_DIR}/others.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(relax python ${CMAKE_CURRENT_SOURCE_DIR}/relax.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(restart python ${CMAKE_CURRENT_SOURCE_DIR}/restart.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(uparams python ${CMAKE_CURRENT_SOURCE_DIR}/uparams.py ${CMAKE_CURRENT_BINARY_DIR})

endif(NOT compile_debug_tests)
