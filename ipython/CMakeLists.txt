set(PYIPY_PY __init__.py _goto.py _explore.py _showme.py _launch.py _collect.py)
install(FILES ${PYIPY_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ipython)

