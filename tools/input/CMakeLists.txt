install( FILES __init__.py block.py keywords.py tree.py
         DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/tools/input )
if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
