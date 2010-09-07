set(PYGAPESCAN_PY __init__.py)
install(FILES ${PYGAPESCAN_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ga/pescan)
add_subdirectory(elemental)

