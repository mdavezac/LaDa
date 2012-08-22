install( FILES __init__.py json.py extract.py makeclass.py uirepr.py block.py keywords.py
         DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/functools )

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
