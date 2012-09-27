set(PYDFTCRYSTAL_PY __init__.py input.py geometry.py basis.py parse.py crystal.py
                    functional.py hamiltonian.py extract.py optgeom.py electronic.py
                    molecule.py relax.py external.py ) 

install(FILES ${PYDFTCRYSTAL_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/dftcrystal)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
