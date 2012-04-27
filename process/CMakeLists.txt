set(PYPROCEED_PY call.py  __init__.py  iterator.py  jobfolder.py 
                 process.py  program.py mpi.py) 
install(FILES ${PYPROCEED_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/process)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
