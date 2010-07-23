set(PYIPY_PY __init__.py _goto.py _explore.py)
install(FILES ${PYIPY_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ipython)

