set(PYGAPESCAN_PY __init__.py functional.py extract.py evaluator.py)
install(FILES ${PYGAPESCAN_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ga/escan
                PERMISSIONS ${PERMISSIONS})
add_subdirectory(elemental)
add_subdirectory(nanowires)

