if(NOT compile_debug_tests)
  set(TYPES  0 1 2 3)
  set(ACCESSES intra external)
  foreach(ACCESS ${ACCESSES})
    foreach(TYPE ${TYPES})
      add_executable(type${ACCESS}${TYPE} types.cc)
      if(${ACCESS} STREQUAL intra)
        set_property( TARGET type${ACCESS}${TYPE} PROPERTY COMPILE_DEFINITIONS
                      LADA_INTERNAL LADA_TEST_TYPE=${TYPE})
      elseif(${ACCESS} STREQUAL external)
        set_property( TARGET type${ACCESS}${TYPE} PROPERTY COMPILE_DEFINITIONS
                      LADA_TEST_STRUCTURE=${STR} LADA_TEST_TYPE=${TYPE})
      endif(${ACCESS} STREQUAL intra)
      ADD_TEST(type${ACCESS}${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/type${ACCESS}${TYPE})
      target_link_libraries(type${ACCESS}${TYPE} lns)
    endforeach(TYPE)
  endforeach(ACCESS)

  add_executable(required_option required_option.cc)
  ADD_TEST(required_option ${CMAKE_CURRENT_BINARY_DIR}/required_option)
  target_link_libraries(required_option lns)

  add_executable(required_section required_section.cc)
  ADD_TEST(required_section ${CMAKE_CURRENT_BINARY_DIR}/required_section)
  target_link_libraries(required_section lns)

  add_executable(id_option id_option.cc)
  ADD_TEST(id_option ${CMAKE_CURRENT_BINARY_DIR}/id_option)
  target_link_libraries(id_option lns)

  add_executable(enum enum.cc)
  ADD_TEST(enum ${CMAKE_CURRENT_BINARY_DIR}/enum)
  target_link_libraries(enum lns)

  add_executable(tie tie.cc)
  ADD_TEST(tie ${CMAKE_CURRENT_BINARY_DIR}/tie)
  target_link_libraries(tie lns)

  set(TYPES  0 1 2)
  foreach(TYPE ${TYPES})
    add_executable(default${TYPE} default.cc)
    target_link_libraries(default${TYPE} lns)
    set_property( TARGET default${TYPE} PROPERTY COMPILE_DEFINITIONS
                  LADA_TEST_TYPE=${TYPE})
    target_link_libraries(default${TYPE} lns)
    ADD_TEST(default${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/default${TYPE})
  endforeach(TYPE)

  set(TYPES  0 1 2)
  foreach(TYPE ${TYPES})
    add_executable(vector_action${TYPE} vector_action.cc)
    set_property( TARGET vector_action${TYPE} PROPERTY COMPILE_DEFINITIONS
                  LADA_TEST_TYPE=${TYPE})
    target_link_libraries(vector_action${ACCESS}${TYPE} lns)
    ADD_TEST(vector_action${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/vector_action${TYPE})
  endforeach(TYPE)

  add_executable(push_back push_back.cc)
  ADD_TEST(push_back ${CMAKE_CURRENT_BINARY_DIR}/push_back)
  target_link_libraries(push_back lns)

  add_executable(merge merge.cc)
  ADD_TEST(merge ${CMAKE_CURRENT_BINARY_DIR}/merge)
  target_link_libraries(merge lns)

  add_executable(shared_ptr shared_ptr.cc)
  ADD_TEST(shared_ptr ${CMAKE_CURRENT_BINARY_DIR}/shared_ptr)
  target_link_libraries(shared_ptr lns)
endif(NOT compile_debug_tests)
