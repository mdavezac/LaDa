set(PYGAPESCAN_PY __init__.py extract.py)
install(FILES ${PYGAPESCAN_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ga/escan)
add_subdirectory(elemental)

