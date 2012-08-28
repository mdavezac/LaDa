install( FILES __init__.py json.py extract.py makeclass.py uirepr.py 
         DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/tools )

add_subdirectory(input)
if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
