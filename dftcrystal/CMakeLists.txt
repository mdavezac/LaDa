set(PYDFTCRYSTAL_PY __init__.py input.py geometry.py)

install(FILES ${PYDFTCRYSTAL_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/dftcrystal)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
