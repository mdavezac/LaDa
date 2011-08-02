set( LIBLNS_HEADERS 
     lns.h string_type.h tags.h access.h parser_base.h
     action/type_to_regex.h action/push_back.h action/fusion.h
     load/load.h 
     save/save.h
     sequencer/binary_const_iterator.h sequencer/binary.h
     sequencer/binary_range.h 
     tree/base.h tree/section.h
     xml/parser.h xml/printer.h xml/tags.h
     xpr/option.h xpr/section.h xpr/section_data.h xpr/subsections.h xpr/parameters.h 
     xpr/merge.h xpr/utilities.h )

set( LIBLNS_SOURCES 
     load/option.cc load/section.cc load/load.cc
     save/save.cc
     sequencer/binary.cc 
     xml/parser.cc xml/printer.cc
     xpr/section.cc )


                  

add_library (lns SHARED ${LIBLNS_SOURCES} ${LIBLNS_HEADERS})
set_target_properties(lns PROPERTIES OUTPUT_NAME lada_lns)
# target_link_libraries(lns opt mathlib)
install(TARGETS lns DESTINATION lib64)

if(NOT noinstall_headers)
  install(FILES ${LIBCRYSTAL_HEADERS} DESTINATION include/LaDa/load_n_save)
endif(NOT noinstall_headers)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
