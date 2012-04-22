add_executable (pifunc functional.cc)
if(NOT compile_debug_tests)
  
  add_executable (pifunc functional.cc)
  add_test(pifunc pifunc --order 2)

endif(NOT compile_debug_tests)

  add_test(program python ${CMAKE_CURRENT_SOURCE_DIR}/program.py ${CMAKE_CURRENT_BINARY_DIR}/pifunc)
  add_test(call python ${CMAKE_CURRENT_SOURCE_DIR}/call.py ${CMAKE_CURRENT_BINARY_DIR}/pifunc)
  add_test(iterator python ${CMAKE_CURRENT_SOURCE_DIR}/iterator.py ${CMAKE_CURRENT_BINARY_DIR}/pifunc)
  add_test(jobfolder python ${CMAKE_CURRENT_SOURCE_DIR}/jobfolder.py ${CMAKE_CURRENT_BINARY_DIR}/pifunc)
