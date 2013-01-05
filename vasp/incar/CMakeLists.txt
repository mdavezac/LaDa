set(PYVASPINCAR_PY __init__.py   _params.py)
install(FILES ${PYVASPINCAR_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/vasp/incar)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
