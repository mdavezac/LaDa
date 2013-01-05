set(PYGAPESCAN_PY __init__.py functional.py converter.py)
install(FILES ${PYGAPESCAN_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/ga/escan/nanowires)

