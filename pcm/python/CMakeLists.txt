set(PYPCM_PY __init__.py)
set(PYPCM_HEADERS  clj.hpp  functional.hpp)

set(PYPCM_SOURCES  clj.cc  functional.cc  python.cc)


include_directories(${PYTHON_INCLUDE_DIRS})
python_add_module (pypcm ${PYPCM_HEADERS} ${PYPCM_SOURCES})
target_link_libraries(pypcm pcm ${Boost_PYTHON_LIBRARY})
mark_as_advanced(PYTHON_ENABLE_MODULE_pypcm)
mark_as_advanced(PYTHON_MODULE_pypcm_BUILD_SHARED)
set_target_properties(pypcm PROPERTIES PREFIX "_")
set_target_properties(pypcm PROPERTIES OUTPUT_NAME pcm)
install(TARGETS pypcm DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/pcm)
install(FILES ${PYPCM_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/pcm)
