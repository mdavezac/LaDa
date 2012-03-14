if(NOT compile_debug_tests)
  
endif(NOT compile_debug_tests)

add_test(common python ${CMAKE_CURRENT_SOURCE_DIR}/common.py ${CMAKE_CURRENT_SOURCE_DIR})
add_test(dft python ${CMAKE_CURRENT_SOURCE_DIR}/dft.py ${CMAKE_CURRENT_SOURCE_DIR})
