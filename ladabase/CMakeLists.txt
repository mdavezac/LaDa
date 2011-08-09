set(LADABASE_PY __init__.py ipython.py extract.py massextract.py misc.py filter.py)
install(FILES ${LADABASE_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ladabase)
