set( LIBLNS_HEADERS 
     string_type.h
     action/enum.h action/string_to_type.h action/type_to_regex.h
     initializer/tagged_options.h initializer/tagged_sections.h
     initializer/section_.h initializer/option_.h
     load/load.h 
     save/save.h
     parser_base/section.h
     sequencer/binary_const_iterator.h sequencer/binary.h
     sequencer/binary_range.h 
     tree/base.h tree/section.h
     xml/parser.h xml/printer.h xml/tags.h
     xpr/option.h xpr/section.h )

set( LIBLNS_SOURCES 
     load/option.cc load/section.cc
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

if(NOT disable_executables)
# add_executable(main_exe main.cc)
# target_link_libraries(main_exe lns)
# add_executable(initializer_exe initializer/main.cc)
  add_executable(load_exe load/main.cc)
  target_link_libraries(load_exe lns)
  add_executable(xml_exe xml/main.cc)
  target_link_libraries(xml_exe lns)
# add_executable(xpr_exe xpr/main.cc)
# target_link_libraries(xpr_exe lns)
endif(NOT disable_executables)
