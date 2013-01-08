if(NOT compile_debug_tests)

  add_library (quantity SHARED "quantity.cc")
  set_target_properties(quantity PROPERTIES PREFIX "_")
  target_link_libraries(quantity ${PYTHON_LIBRARIES})
  add_test(python_quantity python ${CMAKE_CURRENT_SOURCE_DIR}/quantity.py ${CMAKE_CURRENT_BINARY_DIR})

  add_library (pyobject SHARED "pyobject.cc")
  set_target_properties(pyobject PROPERTIES PREFIX "_")
  target_link_libraries(pyobject ${PYTHON_LIBRARIES})
  add_test(python_object python ${CMAKE_CURRENT_SOURCE_DIR}/pyobject.py ${CMAKE_CURRENT_BINARY_DIR})

endif(NOT compile_debug_tests)

