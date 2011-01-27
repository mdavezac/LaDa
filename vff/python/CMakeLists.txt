set(PYVFF_PY __init__.py)
set(PYVFF_HEADERS vff.hpp data.hpp)
set(PYVFF_SOURCES python.cc  vff.cc data.cc)

include_directories(${PYTHON_INCLUDE_DIRS})
python_add_module (pyvff ${PYVFF_HEADERS} ${PYVFF_SOURCES})
target_link_libraries(pyvff vff ${Boost_PYTHON_LIBRARY})
mark_as_advanced(PYTHON_ENABLE_MODULE_pyvff)
mark_as_advanced(PYTHON_MODULE_pyvff_BUILD_SHARED)
set_target_properties(pyvff PROPERTIES PREFIX "_")
set_target_properties(pyvff PROPERTIES OUTPUT_NAME vff)

install(TARGETS pyvff DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/vff)
install(FILES ${PYVFF_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/vff)
