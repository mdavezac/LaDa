set(PYEWALD_PY __init__.py)
set(PYEWALD_F90 ep_com.f90 ewaldf.f90)
set(PYEWALD_HEADERS ewald.h)
set(PYEWALD_SOURCES ewald.cc module.cc erfc.cc)

include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
add_library (pyewald SHARED ${PYEWALD_SOURCES} ${PYEWALD_HEADERS} ${PYEWALD_F90})
set_target_properties(pyewald PROPERTIES OUTPUT_NAME cppwrappers)
set_target_properties(pyewald PROPERTIES PREFIX "")
target_link_libraries(pyewald ${PYTHON_LIBRARY}
                            ${Boost_MATH_C99_LIBRARY} ${Boost_MATH_C99L_LIBRARYl}
                            ${BOOST_C99F_LIBRARY} )
install(TARGETS pyewald DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ewald)
install(FILES ${PYEWALD_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ewald)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
