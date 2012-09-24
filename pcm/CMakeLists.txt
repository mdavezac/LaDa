set(PYPCM_PY __init__.py)
set(PYPCM_F90 ep_com.f90 ewaldf.f90)
set(PYPCM_HEADERS ewald.h)
set(PYPCM_SOURCES ewald.cc module.cc erfc.cc)

include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
add_library (pypcm SHARED ${PYPCM_SOURCES} ${PYPCM_HEADERS} ${PYPCM_F90})
# target_link_libraries(pypcm crystal)
set_target_properties(pypcm PROPERTIES OUTPUT_NAME cppwrappers)
set_target_properties(pypcm PROPERTIES PREFIX "")
target_link_libraries(pypcm ${PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARY} 
                            ${Boost_MATH_C99_LIBRARY} ${Boost_MATH_C99L_LIBRARYl}
                            ${BOOST_C99F_LIBRARY} mathlib)
install(TARGETS pypcm DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/pcm)
install(FILES ${PYPCM_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/pcm)
