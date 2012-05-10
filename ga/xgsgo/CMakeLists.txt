set(PYGA_PY operators.py objective.py random.py __init__.py)
install(FILES ${PYGA_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ga/xgsgo)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
