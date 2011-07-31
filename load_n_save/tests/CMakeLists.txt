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
                    LADA_TEST_STRUCTURE=${STR} LADA_INC_TYPE=${TYPE})
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

add_executable(default default.cc)
ADD_TEST(default ${CMAKE_CURRENT_BINARY_DIR}/default)
target_link_libraries(default lns)



