set(PYIPY_PY __init__.py scattered.py interactive.py asone.py)
install(FILES ${PYIPY_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ipython/launch)

