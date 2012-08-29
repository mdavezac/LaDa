set (PYCE_PY cluster.py __init__.py factory.py)

install(FILES ${PYCE_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ce)
if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
