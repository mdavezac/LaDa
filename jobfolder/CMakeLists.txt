set(PYJOBS_PY __init__.py  runme.py  templates.py runone.py extract.py
              manipulator.py forwarding_dict.py jobfolder.py ordered_dict.py massextract.py)
install(FILES ${PYJOBS_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/jobfolder)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
