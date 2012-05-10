set(PYGA_PY functional.py extract.py standard.py)

configure_file (
 "${PROJECT_SOURCE_DIR}/ga/__init__.in.py"
 "${PROJECT_BINARY_DIR}/ga/__init__.py"
)
install(FILES ${PYGA_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ga)
install(FILES ${PROJECT_BINARY_DIR}/ga/__init__.py DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ga)

add_subdirectory(bitstring)
add_subdirectory(xgsgo)
