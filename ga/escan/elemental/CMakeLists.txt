set(PYGAPESCAN_PY __init__.py evaluator.py extract.py functional.py)
install(FILES ${PYGAPESCAN_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ga/escan/elemental)

