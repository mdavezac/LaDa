set(PYVASP_PY __init__.py extract.py functional.py keywords.py)
install(FILES ${PYVASP_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/gulp)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
