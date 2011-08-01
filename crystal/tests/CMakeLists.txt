# set(TYPES  0 1 2 3)
# set(STRS StructureData TemplateStructure)
# foreach(STR ${STRS})
#   foreach(TYPE ${TYPES})
#     add_executable(init${STR}${TYPE} initialize.cc)
#     set_property( TARGET init${STR}${TYPE} PROPERTY COMPILE_DEFINITIONS
#                   LADA_TEST_STRUCTURE=${STR} LADA_TEST_INCTYPE=${TYPE})
#     ADD_TEST(init${STR}${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/init${STR}${TYPE})
#   endforeach(TYPE)
# endforeach(STR)


# set(TYPES  0 1)
# foreach(TYPE ${TYPES})
#   add_executable(lns_atom${TYPE} lns_atom.cc)
#   set_property( TARGET lns_atom${TYPE} PROPERTY COMPILE_DEFINITIONS
#                 LADA_TEST_INCTYPE=${TYPE})
#   target_link_libraries(lns_atom${TYPE} lns)
#   ADD_TEST(lns_atom${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/lns_atom${TYPE})
# endforeach(TYPE)

set(TYPES  0 1)
set(STRS StructureData TemplateStructure)
foreach(STR ${STRS})
  foreach(TYPE ${TYPES})
    add_executable(lns_structure${STR}${TYPE} lns_structure.cc)
    set_property( TARGET lns_structure${STR}${TYPE} PROPERTY COMPILE_DEFINITIONS
                  LADA_TEST_STRUCTURE=${STR} LADA_TEST_INCTYPE=${TYPE})
    target_link_libraries(lns_structure${STR}${TYPE} lns)
    ADD_TEST(lns_structure${STR}${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/lns_structure${STR}${TYPE})
  endforeach(TYPE)
endforeach(STR)
