set(PYIPY_PY __init__.py scattered.py interactive.py asone.py asone_script.py scattered_script.py)
install(FILES ${PYIPY_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ipython/launch)

