set(PYJOBS_PY __init__.py  runme.py  templates.py runone.py bleeder.py extract.py
              manipulator.py forwarding_dict.py jobdict.py ordered_dict.py massextract.py)
install(FILES ${PYJOBS_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/jobs)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
