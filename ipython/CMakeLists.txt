set(PYIPY_PY __init__.py explore.py savefolders.py goto.py listfolders.py
             showme.py export.py jmol.py copyfolder.py)# qstat.py)
install(FILES ${PYIPY_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ipython)
add_subdirectory(launch)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
