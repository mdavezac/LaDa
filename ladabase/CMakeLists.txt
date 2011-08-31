set(LADABASE_PY __init__.py ipython.py vasp.py misc.py filter.py extracted.py
                groundstates.py massextract.py )
install(FILES ${LADABASE_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ladabase)
