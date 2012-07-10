if(NOT compile_debug_tests)
  add_executable (pifunc functional.cc)
  set_target_properties(pifunc PROPERTIES INCLUDE_DIRECTORIES "${MPI_CXX_INCLUDE_PATH}")
  target_link_libraries(pifunc ${MPI_CXX_LIBRARY})
  
  add_test(process_pifunc ${MPIEXEC} pifunc --order 4)
  add_test(process_program python ${CMAKE_CURRENT_SOURCE_DIR}/program.py ${CMAKE_CURRENT_BINARY_DIR}/pifunc)
  add_test(process_call python ${CMAKE_CURRENT_SOURCE_DIR}/call.py ${CMAKE_CURRENT_BINARY_DIR}/pifunc)
  add_test(process_iterator python ${CMAKE_CURRENT_SOURCE_DIR}/iterator.py ${CMAKE_CURRENT_BINARY_DIR}/pifunc)
  add_test(process_jobfolder python ${CMAKE_CURRENT_SOURCE_DIR}/jobfolder.py ${CMAKE_CURRENT_BINARY_DIR}/pifunc)
  add_test(process_jobfolder python ${CMAKE_CURRENT_SOURCE_DIR}/jobfolder.py ${CMAKE_CURRENT_BINARY_DIR}/pifunc)

endif(NOT compile_debug_tests)
