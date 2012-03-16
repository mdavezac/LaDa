if(NOT compile_debug_tests)
  

endif(NOT compile_debug_tests)
add_test(specie python ${CMAKE_CURRENT_SOURCE_DIR}/species.py ${CMAKE_CURRENT_SOURCE_DIR})
add_test(run python ${CMAKE_CURRENT_SOURCE_DIR}/run.py ${CMAKE_CURRENT_SOURCE_DIR})
add_test(runepi python ${CMAKE_CURRENT_SOURCE_DIR}/runepi.py ${CMAKE_CURRENT_SOURCE_DIR})
add_test(runrelax python ${CMAKE_CURRENT_SOURCE_DIR}/runrelax.py ${CMAKE_CURRENT_SOURCE_DIR})
add_test(runepidoc python ${CMAKE_CURRENT_SOURCE_DIR}/runepidoc.py ${CMAKE_CURRENT_SOURCE_DIR})

