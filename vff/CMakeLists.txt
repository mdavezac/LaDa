set( PYVFF_HEADERS 
     node/pybase.h node/sequence.hpp node/iterator.hpp
     edge/pybase.h )
set( PYVFF_SOURCES 
     module.cc
     node/pybase.cc 
     edge/pybase.cc )
set( PYVFF_PY __init__.py functional.py vff.py _minimize.py )

install(FILES ${PYVFF_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/vff)

include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
add_library (pyvff SHARED ${PYVFF_SOURCES} ${PYVFF_HEADERS})
set_target_properties(pyvff PROPERTIES OUTPUT_NAME cppwrappers)
set_target_properties(pyvff PROPERTIES PREFIX "")
target_link_libraries(pyvff ${PYTHON_LIBRARY} crystal)
install(TARGETS pyvff DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/vff)


if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
