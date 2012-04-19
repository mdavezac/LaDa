set(PYPROCESS_PY call.py  __init__.py  iterator.py  process.py  program.py jobfolder.py)

install(FILES ${PYPROCESS_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/process)

# if(compile_tests OR compile_debug_tests)
#   add_subdirectory(tests)
# endif(compile_tests OR compile_debug_tests)
