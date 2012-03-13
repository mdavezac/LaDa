set(PYVASPEXTRACT_PY __init__.py mixin.py base.py)
install(FILES ${PYVASPEXTRACT_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/vasp/extract)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
