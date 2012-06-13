if(NOT compile_debug_tests)

  add_test(vasp_specie python ${CMAKE_CURRENT_SOURCE_DIR}/species.py ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(vasp_run python ${CMAKE_CURRENT_SOURCE_DIR}/run.py ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(vasp_runepi python ${CMAKE_CURRENT_SOURCE_DIR}/runepi.py ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(vasp_runrelax python ${CMAKE_CURRENT_SOURCE_DIR}/runrelax.py ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(vasp_runepidoc python ${CMAKE_CURRENT_SOURCE_DIR}/runepidoc.py ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(vasp_runmass python ${CMAKE_CURRENT_SOURCE_DIR}/runmass.py ${CMAKE_CURRENT_SOURCE_DIR})

endif(NOT compile_debug_tests)
