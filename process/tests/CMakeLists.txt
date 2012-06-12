if(NOT compile_debug_tests)
  add_executable (pifunc functional.cc)
  set_target_properties(pifunc PROPERTIES INCLUDE_DIRECTORIES ${MPI_CXX_INCLUDE_PATH})
  
  add_test(process_pifunc pifunc --order 2)
  add_test(process_program python ${CMAKE_CURRENT_SOURCE_DIR}/program.py ${CMAKE_CURRENT_BINARY_DIR}/pifunc)
  add_test(process_call python ${CMAKE_CURRENT_SOURCE_DIR}/call.py ${CMAKE_CURRENT_BINARY_DIR}/pifunc)
  add_test(process_iterator python ${CMAKE_CURRENT_SOURCE_DIR}/iterator.py ${CMAKE_CURRENT_BINARY_DIR}/pifunc)
  add_test(process_jobfolder python ${CMAKE_CURRENT_SOURCE_DIR}/jobfolder.py ${CMAKE_CURRENT_BINARY_DIR}/pifunc)
  add_test(process_jobfolder python ${CMAKE_CURRENT_SOURCE_DIR}/jobfolder.py ${CMAKE_CURRENT_BINARY_DIR}/pifunc)

endif(NOT compile_debug_tests)
