set(PYGA_PY __init__.py operators.py)
install(FILES ${PYGA_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ga/bitstring)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
