set(PYGACE_PY __init__.py evaluator.py functional.py)
install(FILES ${PYGACE_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ga/ce)
