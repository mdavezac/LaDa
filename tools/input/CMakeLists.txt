install( FILES __init__.py block.py keywords.py tree.py listblock.py
         DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/tools/input )
if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
