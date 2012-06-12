if(NOT compile_debug_tests)
  
  add_test(incar_algo python ${CMAKE_CURRENT_SOURCE_DIR}/algo.py)
  add_test(incar_choices python ${CMAKE_CURRENT_SOURCE_DIR}/choices.py)
  add_test(incar_magmom python ${CMAKE_CURRENT_SOURCE_DIR}/magmom.py)
  add_test(incar_nelect python ${CMAKE_CURRENT_SOURCE_DIR}/nelect.py)
  add_test(incar_npar python ${CMAKE_CURRENT_SOURCE_DIR}/npar.py)
  add_test(incar_others python ${CMAKE_CURRENT_SOURCE_DIR}/others.py)
  add_test(incar_relax python ${CMAKE_CURRENT_SOURCE_DIR}/relax.py)
  add_test(incar_restart python ${CMAKE_CURRENT_SOURCE_DIR}/restart.py)
  add_test(incar_uparams python ${CMAKE_CURRENT_SOURCE_DIR}/uparams.py)
  add_test(incar_smearing python ${CMAKE_CURRENT_SOURCE_DIR}/uparams.py)

endif(NOT compile_debug_tests)
