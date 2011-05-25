set(PYGAPESCAN_PY __init__.py functional.py extract.py)
install(FILES ${PYGAPESCAN_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ga/escan)
add_subdirectory(elemental)
add_subdirectory(nanowires)

