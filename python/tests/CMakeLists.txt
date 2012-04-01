include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})

if(NOT compile_debug_tests)

  set(TYPES AttributeError KeyError ValueError IndexError TypeError)
  foreach(TYPE ${TYPES})
    python_add_module (exception_${TYPE} "" "exceptions.cc")
    set_property( TARGET exception_${TYPE} PROPERTY COMPILE_DEFINITIONS
                  LADA_TYPE=${TYPE} LADA_TYPENAME="${TYPE}" LADA_MODULE=exception_${TYPE})
    mark_as_advanced(PYTHON_ENABLE_MODULE_exception_${TYPE})
    mark_as_advanced(PYTHON_MODULE_exception_${TYPE}_BUILD_SHARED)
    set_target_properties(exception_${TYPE} PROPERTIES OUTPUT_NAME exception_${TYPE})
    target_link_libraries(exception_${TYPE} ${Boost_PYTHON_LIBRARY})
    configure_file (
     "${CMAKE_CURRENT_SOURCE_DIR}/exceptions_py.py"
     "${CMAKE_CURRENT_BINARY_DIR}/exception_${TYPE}.py"
    )
    add_test(exception_${TYPE} python ${CMAKE_CURRENT_BINARY_DIR}/exception_${TYPE}.py)
  endforeach(TYPE ${TYPES})

  set(TYPES root input internal out_of_range infinite_loop)
  foreach(TYPE ${TYPES})
    python_add_module (exception_${TYPE} "" "exceptions.cc")
    set_property( TARGET exception_${TYPE} PROPERTY COMPILE_DEFINITIONS
                  LADA_TYPE=${TYPE} LADA_TYPENAME="${TYPE}" LADA_MODULE=exception_${TYPE})
    mark_as_advanced(PYTHON_ENABLE_MODULE_exception_${TYPE})
    mark_as_advanced(PYTHON_MODULE_exception_${TYPE}_BUILD_SHARED)
    set_target_properties(exception_${TYPE} PROPERTIES OUTPUT_NAME exception_${TYPE})
    target_link_libraries(exception_${TYPE} ${Boost_PYTHON_LIBRARY})
    configure_file (
     "${CMAKE_CURRENT_SOURCE_DIR}/exceptions.py"
     "${CMAKE_CURRENT_BINARY_DIR}/exception_${TYPE}.py"
    )
    add_test(exception_${TYPE} python ${CMAKE_CURRENT_BINARY_DIR}/exception_${TYPE}.py)
  endforeach(TYPE ${TYPES})

  add_library (quantity SHARED "quantity.cc")
  target_link_libraries(quantity ${PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARY})

  set_target_properties(quantity PROPERTIES PREFIX "_")
  add_test(quantity python ${CMAKE_CURRENT_SOURCE_DIR}/quantity.py ${CMAKE_CURRENT_BINARY_DIR})

endif(NOT compile_debug_tests)
