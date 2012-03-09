if(NOT compile_debug_tests)
endif(NOT compile_debug_tests)
  
  add_test(algo python ${CMAKE_CURRENT_SOURCE_DIR}/algo.py)
  add_test(choices python ${CMAKE_CURRENT_SOURCE_DIR}/choices.py)
  add_test(magmom python ${CMAKE_CURRENT_SOURCE_DIR}/magmom.py)
  add_test(nelect python ${CMAKE_CURRENT_SOURCE_DIR}/nelect.py)
  add_test(npar python ${CMAKE_CURRENT_SOURCE_DIR}/npar.py)
  add_test(others python ${CMAKE_CURRENT_SOURCE_DIR}/others.py)
  add_test(relax python ${CMAKE_CURRENT_SOURCE_DIR}/relax.py)
  add_test(restart python ${CMAKE_CURRENT_SOURCE_DIR}/restart.py)
  add_test(uparams python ${CMAKE_CURRENT_SOURCE_DIR}/uparams.py)
  add_test(smearing python ${CMAKE_CURRENT_SOURCE_DIR}/uparams.py)

