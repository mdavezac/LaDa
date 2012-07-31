set(PYVASP_PY files.py __init__.py specie.py functional.py relax.py emass.py
              keyword.py )

install(FILES ${PYVASP_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/vasp)

# add_subdirectory(incar)
add_subdirectory(extract)
add_subdirectory(nlep)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
